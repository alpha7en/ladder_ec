"""Instance checks (Sec. 7.1, Table 5): fixed-size checks of the kind used to grade such tasks today.

judge_state ("State"): n = 1..10; bits a: all for n <= 3, otherwise 0, 1...1 and 3 random strings; compares the
  outputs on |0...0>.
judge_unitary ("Rand."): the same n and a; compares the action on 2 random product states.
random_test ("Samp."): 20 cases, n uniform in 1..20, a uniform in [0, 2^n), one random product state per case;
  cases where the reference fails are skipped.
Comparison by the task relation p["label"], tolerance 1e-9:
  E2: |<r|g>| = 1 for circuits of equal width (different widths: caught);
  E5phase: the same, but all three checks use only the input |0...0>, the relation speaks only about it;
  E3: ancillas in the most significant qubits (circuit width minus n), input data (x) |0>; the solution must return
  the ancillas to |0> and agree with the reference on the data exactly; ancilla counts may differ.
A solution failing where the reference runs is caught. A case with a circuit wider than WMAX qubits is not
simulated and is skipped (field "wide"). A check stops at the first caught case, as a judge does.
State vectors only, no full matrices. Random generator of a pair: numpy.random.default_rng(pair number); from it, in
this order, the bit strings of the judges (shared by both), the plan of random_test, then the input states. Pairs
with a no-op mutation and reference-reference pairs are not checked. One process (the times of Sec. 7.4 are
measured on the machine of the method); an interrupted run resumes.
Usage: python -I experiments/instance_checks.py   -> data/runs/instance_checks.jsonl
"""
import json
import multiprocessing
import os
import sys
import time
import numpy as np
from qiskit.quantum_info import Statevector

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "experiments", "benchmark"))
from ground_truth import build, fn, product_state   # noqa: E402

NJ, NR, NR_MAX = 10, 20, 20
WMAX = 20          # 2^20 amplitudes, 16 MB; random_test goes up to n = 20 and must not skip circuits of width n


def zero(w):
    return np.eye(1, 2 ** w)[0].astype(complex)


def out(C, v, w):
    """output of circuit C on input v (on w qubits) padded with ancillas |0> in the most significant qubits"""
    return Statevector(np.kron(zero(C.num_qubits - w), v)).evolve(C).data


def differs(label, n, R, G, inputs):
    """True: the solution differs from the reference on at least one input (G = None: the solution failed);
    None: the case is wider than WMAX. inputs(w): inputs on w qubits"""
    if G is None:
        return True
    wR, wG = R.num_qubits, G.num_qubits
    if label != "E3" and wR != wG:
        return True
    if max(wR, wG) > WMAX:
        return None
    if label != "E3":
        return any(abs(abs(np.vdot(out(R, v, wR), out(G, v, wR))) - 1) >= 1e-9 for v in inputs(wR))
    d = 2 ** n
    for v in inputs(n):
        r, g = out(R, v, n)[:d], out(G, v, n)[:d]
        if abs(np.linalg.norm(g) - 1) >= 1e-9 or not np.allclose(r, g, atol=1e-9):
            return True
    return False


def run(label, cases, inputs, fR, fG):
    """cases: (n, a) in order; inputs(w): input states. Returns the record of one instance check."""
    t0, k, skipped, wide = time.perf_counter(), 0, 0, 0
    first = None
    for n, a in cases:
        R = build(fR, n, a)
        if R is None:
            skipped += 1
            continue
        d = differs(label, n, R, build(fG, n, a), inputs)
        if d is None:
            wide += 1
            continue
        k += 1
        if d:
            first = dict(n=n, a=a)
            break
    return dict(caught=first is not None, first=first, instances=k, skipped=skipped, wide=wide,
                time_s=round(time.perf_counter() - t0, 4))


def check(job):
    t, p = job
    rng = np.random.default_rng(t)
    label = p.get("label", "E2")
    fR, fG = fn(p["ref"]), fn(p["sol"])
    bits = fR.__code__.co_argcount >= 2
    judge = [(n, a) for n in range(1, NJ + 1) for a in (
        (range(2 ** n) if n <= 3 else sorted({0, 2 ** n - 1} | {int(x) for x in rng.integers(0, 2 ** n, 3)}))
        if bits else [None])]
    plan = []
    for _ in range(NR):
        n = int(rng.integers(1, NR_MAX + 1))
        plan.append((n, int(rng.integers(0, 2 ** n)) if bits else None))
    state = lambda w: [zero(w)]                                                  # noqa: E731
    prod = lambda k: (lambda w: [np.asarray(product_state(w, rng)) for _ in range(k)])  # noqa: E731
    r = dict(id=p["id"], label=label)
    r["judge_state"] = run(label, judge, state, fR, fG)
    r["judge_unitary"] = run(label, judge, state if label == "E5phase" else prod(2), fR, fG)
    r["random_test"] = run(label, plan, state if label == "E5phase" else prod(1), fR, fG)
    return r


def main():
    bench = os.path.join(ROOT, "data", "benchmark")
    pairs = [json.loads(l) for l in open(os.path.join(bench, "pairs.jsonl"))]
    truth = {r["id"]: r for r in map(json.loads, open(os.path.join(bench, "truth.jsonl")))}
    out = os.path.join(ROOT, "data", "runs", "instance_checks.jsonl")
    done = {json.loads(l)["id"] for l in open(out)} if os.path.exists(out) else set()
    # reference-reference pairs (style:self, tasks outside the class) are skipped: nothing to catch, and mcx of
    # large arity up to n = 20 takes hours to simulate
    jobs = [(t, p) for t, p in enumerate(pairs) if not truth[p["id"]]["noop"] and p["id"] not in done
            and p["kind"] != "style:self"]
    t0 = time.time()
    with open(out, "a") as f, multiprocessing.Pool(1) as pool:
        for k, r in enumerate(pool.imap_unordered(check, jobs)):
            f.write(json.dumps(r, ensure_ascii=False) + "\n")
            f.flush()
            if k % 50 == 0:
                print(k, "of", len(jobs), "%.0f s" % (time.time() - t0), flush=True)
    print("run %.0f s" % (time.time() - t0), flush=True)


if __name__ == "__main__":
    main()
