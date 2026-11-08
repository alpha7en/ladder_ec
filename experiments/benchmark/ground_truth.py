"""Ground truth of the benchmark by direct computation with Qiskit (Sec. 7.1).

Without bits: n = 1..7, full matrices (qiskit Operator). With bits: n = 1..6, all a for n <= 4, otherwise 0, 1...1
and 6 random strings. In addition the sizes check_n of each pair (computed by make_benchmark.py: c + 2 and c + 3 for
loop starts c < 10, K + 1 and K + 2 for thresholds K, and the sizes where a loop changes its number of iterations,
up to n = 24), there by the state vector
on |0...0> and on a random product state; with bits at a = 0, 1...1, 0101..., 1010... and one random string.
Relations: E1, E2 (up to a phase, its own for every n, a), E5, E5 up to phase; for tasks with relation E3 only E3
and E3 up to phase (ancillas: the last qubits beyond n, |0> on input, both circuits must return them to |0>).
Circuits wider than 10 qubits (tasks outside the class) are checked by the state vector only. A solution that
fails where the reference runs differs from it at that n. Sizes go upward; once every relation of a pair is
"not equal" the search stops (the truth and the first difference do not change; the row then has the field
"stopped" with that size). One row per pair in data/benchmark/truth.jsonl: truth, first_difference and noop
(the edit did not change the circuit at any checked n; such pairs are left out of the tables). An interrupted
run resumes; three processes.
Usage: python -I experiments/benchmark/ground_truth.py
"""
import importlib.util
import json
import multiprocessing
import os
import numpy as np
from qiskit.quantum_info import Operator, Statevector

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
BENCH = os.path.join(ROOT, "data", "benchmark")
NB, NB_BITS, WIDE = 7, 6, 10
_mods = {}


def fn(ref):
    path = os.path.join(ROOT, ref["file"])
    if path not in _mods:
        spec = importlib.util.spec_from_file_location(os.path.basename(path)[:-3], path)
        _mods[path] = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(_mods[path])
    return getattr(_mods[path], ref["function"])


def build(f, n, a):
    try:
        return f(n, a) if a is not None else f(n)
    except Exception:
        return None


_refs = {}
NONE = dict(E1=False, E2=False, E5=False, E5phase=False)


def clean(v, d):
    """the ancillas returned to |0>: zeros outside the first d amplitudes"""
    return np.allclose(v[d:], 0, atol=1e-9)


def compare_full(key, R, G, n):
    if G is None:
        return NONE
    if key not in _refs:
        _refs[key] = Operator(R).data
    UR, UG = _refs[key], Operator(G).data
    if UR.shape != UG.shape:
        return NONE
    d = len(UR)
    return dict(E1=np.allclose(UR, UG, atol=1e-9), E2=abs(abs(np.vdot(UR, UG)) - d) < 1e-7,
                E5=np.allclose(UR[:, 0], UG[:, 0], atol=1e-9), E5phase=abs(abs(np.vdot(UR[:, 0], UG[:, 0])) - 1) < 1e-9)


def compare_e3(key, R, G, n):
    """E3: G(psi (x) 0) = (V psi) (x) 0 with V the part of the reference on the data; ancilla counts may differ"""
    if G is None:
        return dict(E3=False, E3phase=False)
    if key not in _refs:
        _refs[key] = Operator(R).data
    d, UR, UG = 2 ** n, _refs[key], Operator(G).data
    assert all(clean(UR[:, x], d) for x in range(d)), "reference does not return its ancillas"
    ok = all(clean(UG[:, x], d) for x in range(d))
    V, W = UR[:d, :d], UG[:d, :d]
    return dict(E3=ok and np.allclose(V, W, atol=1e-9), E3phase=ok and abs(abs(np.vdot(V, W)) - d) < 1e-7)


def product_state(n, rng):
    v = np.ones(1)
    for _ in range(n):
        q = rng.normal(size=2) + 1j * rng.normal(size=2)
        v = np.kron(q / np.linalg.norm(q), v)
    return v


def padded(v, qc):
    """data state with the ancillas of the circuit (qubits beyond len(v)) in |0>"""
    e = np.zeros(2 ** qc.num_qubits // len(v))
    e[0] = 1
    return Statevector(np.kron(e, v))


def compare_sv(R, G, n, rng, e3):
    if G is None:
        return dict(E3=False, E3phase=False) if e3 else NONE
    if e3:
        d, out = 2 ** n, dict(E3=True, E3phase=True)
        for v in (np.eye(1, d)[0], product_state(n, rng)):
            r, g = padded(v, R).evolve(R).data, padded(v, G).evolve(G).data
            assert clean(r, d), "reference does not return its ancillas"
            ok = clean(g, d)
            out["E3"] &= ok and np.allclose(r[:d], g[:d], atol=1e-9)
            out["E3phase"] &= ok and abs(abs(np.vdot(r[:d], g[:d])) - 1) < 1e-9
        return out
    if R.num_qubits != G.num_qubits:
        return NONE
    w = R.num_qubits
    z = Statevector.from_int(0, 2 ** w)
    r0, g0 = z.evolve(R).data, z.evolve(G).data
    psi = Statevector(product_state(w, rng))
    r, g = psi.evolve(R).data, psi.evolve(G).data
    e2 = abs(abs(np.vdot(r, g)) - 1) < 1e-9
    return dict(E1=e2 and np.allclose(r, g, atol=1e-9), E2=e2, E5=np.allclose(r0, g0, atol=1e-9),
                E5phase=abs(abs(np.vdot(r0, g0)) - 1) < 1e-9)


def strings(n, rng):
    if n <= 4:
        return range(2 ** n)
    return sorted({0, 2 ** n - 1} | {int(x) for x in rng.integers(0, 2 ** n, 6)})


def check(job):
    t, p = job
    rng = np.random.default_rng(t)
    fR, fG = fn(p["ref"]), fn(p["sol"])
    bits = fR.__code__.co_argcount >= 2
    e3 = p["label"].startswith("E3")
    truth = dict(E3=True, E3phase=True) if e3 else dict(E1=True, E2=True, E5=True, E5phase=True)
    first = {}

    def note(res, where):
        for k, v in res.items():
            if not v and truth[k]:
                truth[k], first[k] = False, where

    sizes = range(1, (NB_BITS if bits else NB) + 1)
    stopped = None
    for n in sizes:
        if not any(truth.values()):
            stopped = n
            break
        for a in (strings(n, rng) if bits else [None]):
            R = build(fR, n, a)
            if R is None:
                continue
            G = build(fG, n, a)
            if R.num_qubits > WIDE:
                note(compare_sv(R, G, n, rng, e3), dict(n=n, a=a, statevector=True))
            else:
                note((compare_e3 if e3 else compare_full)((p["ref"]["function"], n, a), R, G, n), dict(n=n, a=a))
    for n in p["check_n"]:
        if not any(truth.values()):
            stopped = stopped or n
            break
        alt = int("01" * n, 2) & (2 ** n - 1)                          # 0101...: neighbouring bits differ
        for a in ([0, 2 ** n - 1, alt, alt ^ (2 ** n - 1), int(rng.integers(0, 2 ** n))] if bits else [None]):
            R = build(fR, n, a)
            note(compare_sv(R, build(fG, n, a), n, rng, e3), dict(n=n, a=a, statevector=True))
    truth.update({k: None for k in ("E1", "E2", "E5", "E5phase", "E3", "E3phase") if k not in truth})
    truth["checked"] = "n <= %d%s%s" % (
        sizes[-1], ", a: all for n <= 4, 8 for n > 4" if bits else "", "".join(", statevector n = %d" % n for n in p["check_n"]))
    L = p["label"]
    noop = p["expect"][L] is False and truth[L]
    if p["expect"][L] is True and not truth[L]:
        print("MISMATCH", p["id"], first, flush=True)
    row = dict(id=p["id"], truth=truth, first_difference=first, noop=noop)
    if stopped:
        row["stopped"] = stopped
    return row


def main():
    pairs = [json.loads(l) for l in open(os.path.join(BENCH, "pairs.jsonl"))]
    out = os.path.join(BENCH, "truth.jsonl")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    done = {json.loads(l)["id"] for l in open(out)} if os.path.exists(out) else set()
    jobs = [(t, p) for t, p in enumerate(pairs) if p["id"] not in done]
    with open(out, "a") as f, multiprocessing.Pool(3) as pool:
        for k, r in enumerate(pool.imap_unordered(check, jobs)):
            f.write(json.dumps(r, ensure_ascii=False) + "\n")
            f.flush()
            if k % 100 == 0:
                print(k, "of", len(jobs), flush=True)
    rows = [json.loads(l) for l in open(out)]
    print(len(rows), "pairs,", sum(r["noop"] for r in rows), "noop")


if __name__ == "__main__":
    main()
