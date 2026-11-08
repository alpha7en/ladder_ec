"""QCEC enumeration of the sizes n_min <= n <= cutoff (Corollary 1): for tasks with bits all a (if 2^n <= 4096,
otherwise 4096 random strings, marked "a_sampled"), a limit per instance and per pair.

python -I experiments/qcec/driver.py --mode A --pairs PAIRS.jsonl [--pairs ...] [--only ID ...] [--jobs K]
       [--out results.jsonl] [--inst instances.jsonl] [--resume]

Mode A is the enumeration used in the paper; mode B (n = 2..8, 16 random a) imitates a judge and is not reported.
The enumeration of a pair stops at the first "not equal" (a counterexample). An instance that timed out prevents the
answer "equal" for the pair: the result is then "unknown"; a pair over its limit is "pair_timeout"; a pair whose
every instance failed to build is "error".
The cutoff comes from the pair row: field "cutoff" (a number, or a dict relation -> number); n_min from the field
"n_min" (default 1). Each pair runs in its own worker process (worker.py); one line per instance goes to --inst,
the result per pair to --out.
"""
import argparse
import hashlib
import inspect
import json
import os
import random
import resource
import select
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
PY = os.environ.get("QCEC_PY", sys.executable)
INST_TIMEOUT = 60.0
PAIR_LIMIT = 2 * 3600.0
MAX_A = 4096
JUDGE_A = 16
JUDGE_N = range(2, 9)
MEM_LIMIT = 6 << 30
_lock = threading.Lock()

CONFIG = {
    "tool": "mqt.qcec 3.11.0 (mqt-core 3.11.0), qiskit 2.5.2",
    "E1_E2": "verify(run_zx_checker=False): alternating + simulation, default tolerance",
    "E5": "all qubits are ancillas in |0>, without check_partial_equivalence (that one compares probabilities), ZX off",
    "E3": "padding with qubits, qubits >= min(number of qubits) are ancillas in |0>; E3 exact, E3phase up to global phase",
    "timeout_instance_s": INST_TIMEOUT, "limit_pair_s": PAIR_LIMIT,
    "transpile": "basis BASIS of experiments/qcec/worker.py, optimization_level=0, qubits_initially_zero=False; global_phase moved into the gates P(t) X P(t) X",
}


def seed_for(*parts):
    return int(hashlib.sha256("|".join(map(str, parts)).encode()).hexdigest()[:12], 16)


def has_bits(ref):
    """whether the generator takes a second argument (bits a); read from the signature, without running qiskit"""
    import ast
    path = ref["file"] if os.path.isabs(ref["file"]) else os.path.join(ROOT, ref["file"])
    tree = ast.parse(open(path).read())
    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef) and node.name == ref["function"]:
            return len(node.args.args) >= 2
    raise KeyError(ref["function"])


def base_label(pair):
    lab = pair.get("label", "E2")
    lab = lab[0] if isinstance(lab, list) else lab
    return lab.split()[0]


def cutoffs(pair):
    """relation -> cutoff in n, from the field cutoff (a number or a dict)"""
    c = pair.get("cutoff")
    if isinstance(c, dict):
        return {k: int(v) for k, v in c.items() if v is not None}
    if c is None:
        return {}
    return {base_label(pair): int(c)}


def has_bits_pair(pair):
    return has_bits(pair["ref"]) or has_bits(pair["sol"])


def e3_free(pair):
    """number of free input qubits for E3 when the task has a restricted input (field "e3_free" of the pair);
    None: for E3 the ancillas are the qubits beyond the smaller qubit count"""
    return pair.get("e3_free")


def a_values(n, bits, mode, key):
    if not bits:
        return [None], False
    total = 1 << n
    k = MAX_A if mode == "A" else JUDGE_A
    if total <= k:
        return list(range(total)), False
    rng = random.Random(seed_for(key, n))
    return sorted(rng.sample(range(total), k)), True


class Worker:
    def __init__(self):
        self.start()

    def start(self):
        def lim():
            resource.setrlimit(resource.RLIMIT_AS, (MEM_LIMIT, MEM_LIMIT))
        self.p = subprocess.Popen([PY, "-I", os.path.join(HERE, "worker.py")], stdin=subprocess.PIPE,
                                  stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True, preexec_fn=lim)
        self.readline(300)

    def readline(self, timeout):
        r, _, _ = select.select([self.p.stdout], [], [], timeout)
        if not r:
            return None
        line = self.p.stdout.readline()
        return json.loads(line) if line else {"error": "worker process exited (memory?)"}

    def ask(self, job, timeout):
        self.p.stdin.write(json.dumps(job) + "\n")
        self.p.stdin.flush()
        out = self.readline(timeout)
        if out is None or "error" in out and "exited" in out.get("error", ""):
            self.kill()
            self.start()
            return out or {"verdict": "killed", "error": f"no answer within {timeout} s"}
        return out

    def kill(self):
        try:
            self.p.kill()
            self.p.wait(5)
        except Exception:
            pass


def run_pair(pair, crit, cut, mode, inst_path):
    key = f"{pair['id']}|{crit}|{mode}"
    n_min = pair.get("n_min") or 1
    n_min = int(n_min)
    bits = has_bits_pair(pair)
    free = e3_free(pair) if crit in ("E3", "E3phase") else None
    # sizes n where the generator cannot build the circuit are recorded in n_undefined
    ns = list(range(n_min, cut + 1)) if mode == "A" else list(JUDGE_N)
    w = Worker()
    t_start = time.perf_counter()
    summary = {"id": pair["id"], "task": pair.get("task"), "origin": pair.get("origin"),
               "mutation": pair.get("mutation"), "mode": mode, "crit": crit, "cutoff": cut, "n_min": n_min,
               "bits": bits, "e3_free": free, "n_planned": [ns[0], ns[-1]] if ns else None}
    inst = timeouts = errors = 0
    sampled = False
    sum_build = sum_check = 0.0
    verdict, cex, n_done, prefix_ok = "equal", None, None, True
    any_ok, undefined, build_err = False, set(), 0
    try:
        for n in ns:
            avals, smp = a_values(n, bits, mode, key)
            sampled |= smp
            n_ok = True
            for a in avals:
                if time.perf_counter() - t_start > PAIR_LIMIT:
                    verdict = "pair_timeout"
                    break
                job = {"ref": pair["ref"], "sol": pair["sol"], "n": n, "a": a, "crit": crit, "timeout": INST_TIMEOUT, "e3_free": free}
                out = w.ask(job, INST_TIMEOUT * 3 + 120)
                inst += 1
                rec = {"id": pair["id"], "mode": mode, "crit": crit, "n": n, "a": a, **out}
                with _lock, open(inst_path, "a") as f:
                    f.write(json.dumps(rec, ensure_ascii=False) + "\n")
                sum_build += out.get("build_s", 0) or 0
                sum_check += out.get("check_s", 0) or 0
                v = out.get("verdict")
                if "build_error" in out:
                    build_err += 1          # the generator cannot build the circuit at (n, a): not an answer of QCEC
                    if not any_ok:
                        undefined.add(n)
                    continue
                if v in ("equal", "not_equal"):
                    any_ok = True
                if v == "not_equal":
                    verdict, cex = "not_equal", {"n": n, "a": a, "qcec_result": out.get("qcec_result")}
                    break
                if v != "equal":
                    n_ok = False
                    if "error" in out and v is None:
                        errors += 1
                    else:
                        timeouts += 1
            if verdict in ("not_equal", "pair_timeout"):
                break
            if n in undefined:
                continue
            if n_ok and prefix_ok:
                n_done = n
            else:
                prefix_ok = False
    finally:
        w.kill()
    if verdict == "equal" and (timeouts or errors):
        verdict = "unknown"
    if not any_ok and verdict == "equal":
        verdict = "error"           # no instance checked (all failed to build)
    summary.update({"verdict": verdict, "counterexample": cex, "n_all_equal_upto": n_done,
                    "instances": inst, "n_undefined": sorted(undefined), "timeouts": timeouts, "errors": errors, "build_errors": build_err, "a_sampled": sampled,
                    "build_s": round(sum_build, 3), "check_s": round(sum_check, 3),
                    "wall_s": round(time.perf_counter() - t_start, 3), "config": CONFIG})
    return summary


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["A", "B"], required=True)
    ap.add_argument("--pairs", action="append", required=True)
    ap.add_argument("--only", nargs="*")
    ap.add_argument("--jobs", type=int, default=2)
    ap.add_argument("--out", default="results_qcec.jsonl")
    ap.add_argument("--inst", default="results_qcec_instances.jsonl")
    ap.add_argument("--resume", action="store_true", help="skip (pair, relation, mode) already in --out")
    args = ap.parse_args()
    CONFIG["parallel_pairs"] = args.jobs
    CONFIG["machine"] = f"{os.cpu_count()} cores, cloud container"
    pairs = [json.loads(l) for p in args.pairs for l in open(p) if l.strip()]
    if args.only:
        pairs = [p for p in pairs if p["id"] in args.only]
    tasks = []
    done = set()
    if args.resume and os.path.exists(args.out):
        for l in open(args.out):
            r = json.loads(l)
            done.add((r["id"], r.get("crit"), r["mode"]))
    for p in pairs:
        cs = cutoffs(p)
        if not cs and args.mode == "A":
            if (p["id"], None, "A") in done:
                continue
            with _lock, open(args.out, "a") as f:
                f.write(json.dumps({"id": p["id"], "task": p.get("task"), "origin": p.get("origin"), "mode": "A", "verdict": "skipped",
                                    "reason": "no cutoff"}, ensure_ascii=False, default=str) + "\n")
            continue
        if not cs:
            cs = {base_label(p): None}
        for crit, cut in cs.items():
            tasks.append((p, crit, cut))
    if args.resume:
        tasks = [t for t in tasks if (t[0]["id"], t[1], args.mode) not in done]

    def go(t):
        p, crit, cut = t
        try:
            s = run_pair(p, crit, cut, args.mode, args.inst)
        except Exception as e:
            s = {"id": p["id"], "mode": args.mode, "crit": crit, "verdict": "error", "error": repr(e)[:300]}
        with _lock, open(args.out, "a") as f:
            f.write(json.dumps(s, ensure_ascii=False) + "\n")
        print(s["id"], crit, s.get("verdict"), s.get("wall_s"), flush=True)

    with ThreadPoolExecutor(args.jobs) as ex:
        list(ex.map(go, tasks))


if __name__ == "__main__":
    main()
