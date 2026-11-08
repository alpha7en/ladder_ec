"""Exact plan without the cap T <= T_max (Sec. 7.4, Fig. S1): for every pair of the main run answered "equal" with
an error bound (the probabilistic plan, Thm. S.13), the probabilistic plan as in the run (T_max = 8) and the exact
plan of Thm. S.12 with no cap on T, both timed on the automata part only (decide). The exact plan certifies "equal"
with error 0. Primes and n_cut of E2 use the bound of Table 1 for the product series, as decide does.
Input data/runs/method_main.jsonl; output data/runs/exact_plan.jsonl (resumable, one line per pair, written and
synced at once). Limit per pair and plan: 1800 s by default.
Usage: python -I experiments/exact_plan.py [limit_s]"""
import json, os, signal, sys, time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from ladder_ec.check import KINDS, family      # noqa: E402
from ladder_ec.translator import analyze        # noqa: E402
from ladder_ec.plans import plan                # noqa: E402
from ladder_ec.field import field_order         # noqa: E402
from ladder_ec.decide import angles, decide     # noqa: E402

NO_CAP = 10 ** 9


class Timeout(Exception):
    pass


def timed(f, limit):
    signal.alarm(limit)
    t0 = time.time()
    try:
        return f(), round(time.time() - t0, 3)
    except Timeout:
        return None, round(time.time() - t0, 3)
    finally:
        signal.alarm(0)


def run(p, limit):
    ref, sol = [(os.path.join(ROOT, p[s]["file"]), p[s]["function"]) for s in ("ref", "sol")]
    G, R = family(analyze(*sol)), family(analyze(*ref))
    kind = KINDS[p["label"]]
    m = field_order(angles(G) + angles(R))
    row = {"id": p["id"], "label": p["label"], "m": m}
    out, row["random_s"] = timed(lambda: decide(G, R, kind), limit)
    row["random"] = out and dict(verdict=out[0], primes=out[1].get("primes"), error=out[1].get("error"),
                                 n_cut=out[1].get("n_cut"))
    if kind in ("E3", "E3p"):                 # the plan of decide_e3: words on the data qubits
        args = ((G, R), m, kind == "E3p", max(G.n_min, R.n_min) + max(G.extra, R.extra), 1)
    else:
        args = ((G, R), m, "product" if kind == "E2" else kind == "E5p", max(G.n_min, R.n_min) + G.extra, 1)
    pl, row["plan_s"] = timed(lambda: plan(*args, 40, 0, NO_CAP), limit)
    row["T"] = pl and len(pl[0])
    row["plan_exact"] = pl and pl[1] == 0           # False: T would exceed NO_CAP, the plan is probabilistic
    out, row["exact_s"] = timed(lambda: decide(G, R, kind, exact_max=NO_CAP), limit) if pl else (None, None)
    row["exact"] = out and dict(verdict=out[0], primes=out[1].get("primes"), error=out[1].get("error"),
                                n_cut=out[1].get("n_cut"))
    return row


def main():
    limit = int(sys.argv[1]) if len(sys.argv) > 1 else 1800
    signal.signal(signal.SIGALRM, lambda *_: (_ for _ in ()).throw(Timeout()))
    pairs = {p["id"]: p for p in map(json.loads, open(os.path.join(ROOT, "data", "benchmark", "pairs.jsonl")))}
    rows = [json.loads(l) for l in open(os.path.join(ROOT, "data", "runs", "method_main.jsonl"))]
    todo = [pairs[r["id"]] for r in rows if isinstance(r["ours"].get("verdicts"), dict)
            and any("error" in v for v in r["ours"]["verdicts"].values())]
    out = os.path.join(ROOT, "data", "runs", "exact_plan.jsonl")
    done = {json.loads(l)["id"] for l in open(out)} if os.path.exists(out) else set()
    for p in [p for p in todo if p["id"] not in done]:
        row = run(p, limit)
        with open(out, "a") as f:          # written and fsynced line by line
            f.write(json.dumps(row) + "\n")
            f.flush()
            os.fsync(f.fileno())
        print(row["id"], row["T"], row["random_s"], row["exact_s"], row["exact"] and row["exact"]["verdict"], flush=True)


if __name__ == "__main__":
    main()
