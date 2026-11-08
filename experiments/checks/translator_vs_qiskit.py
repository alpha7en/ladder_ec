"""Gate-by-gate check of the translator against Qiskit (Sec. 6.1).

The translator is static: it never runs the generator. This check runs it. For every function of the benchmark
that the translator accepts, the ladder program is expanded (translator.expand) and compared, gate by gate, with the
circuit Qiskit builds, at five sizes n >= n_min and four random bit strings per size.
The functions circuit_gates, same and diag are also used by experiments/external/quantumkatas/run_katas.py.

Usage: python -I experiments/checks/translator_vs_qiskit.py      (prints the counts and the first mismatches)
"""
import collections
import importlib.util
import json
import math
import os
import random
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, ROOT)
from ladder_ec.translator import analyze, expand   # noqa: E402


def circuit_gates(qc):
    """(name, qubit indices, first parameter or None) for every instruction of a Qiskit circuit"""
    out = []
    for inst in qc.data:
        op = inst.operation
        out.append((op.name, tuple(qc.find_bit(q).index for q in inst.qubits),
                    float(op.params[0]) if op.params else None))
    return out


def same(a, b):
    if len(a) != len(b):
        return False
    for (o1, q1, a1), (o2, q2, a2) in zip(a, b):
        if o1 != o2 or q1 != q2:
            return False
        if (a1 is None) != (a2 is None) or (a1 is not None and not math.isclose(a1, a2, abs_tol=1e-12)):
            return False
    return True


def access_kinds(desc):
    """how each bit parameter is read: "index" (a list, a[k]) or "shift" (an integer, (a >> k) & 1)"""
    kinds = {}
    for seg in desc["segments"]:
        for g in seg.get("gates") or seg.get("body"):
            for c in g.get("cond", []):
                kinds[c["param"]] = c["access"]
    return kinds


def diag(desc, fn, n_extra=5, samples=4, seed=0):
    """compare the expansion of the ladder program desc with fn(n, ...) at n_min <= n < n_min + n_extra;
    bit strings have length n + 4. Returns (True, "<k> cases") or (False, "<first failing n and a>")."""
    rng = random.Random(seed)
    kinds = access_kinds(desc)
    checked = 0
    for n in range(desc["n_min"], desc["n_min"] + n_extra):
        for _ in range(samples if desc["bit_params"] else 1):
            params = {}
            for p in desc["bit_params"]:
                width = n + 4
                params[p] = [rng.randint(0, 1) for _ in range(width)] if kinds.get(p) == "index" else \
                    rng.randrange(2 ** width)
            if not same(circuit_gates(fn(n, **params)), expand(desc, n, params)):
                return False, f"n={n}, params={params}"
            checked += 1
    return True, f"{checked} cases"


def main():
    funcs = {}
    for line in open(os.path.join(ROOT, "data", "benchmark", "pairs.jsonl")):
        p = json.loads(line)
        for side in ("ref", "sol"):
            funcs[(p[side]["file"], p[side]["function"])] = 1
    mods, st, bad = {}, collections.Counter(), []
    for file, fn in sorted(funcs):
        path = os.path.join(ROOT, file)
        d = analyze(path, fn)
        if d["status"] != "ladder":
            st["refused"] += 1
            continue
        if file not in mods:
            spec = importlib.util.spec_from_file_location(os.path.basename(path)[:-3], path)
            mods[file] = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(mods[file])
        ok, info = diag(d, getattr(mods[file], fn))
        st["ok" if ok else "mismatch"] += 1
        if not ok:
            bad.append((file, fn, info))
    print(len(funcs), "functions", dict(st))
    for b in bad[:30]:
        print("MISMATCH", *b)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
