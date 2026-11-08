"""Shared part of ladder_dims.py and theorem1_bound.py: the distinct loops of the benchmark ladders and the shape
parameters of each loop that enter the bound of Theorem 1 (Sec. 5.2, Table 1)."""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, ROOT)
from ladder_ec.translator import analyze     # noqa: E402
from ladder_ec.check import family           # noqa: E402
from ladder_ec.ladders import Loop, mirror   # noqa: E402


def loops():
    """yield (function, file, family, loop segment, shape) for every distinct loop of the accepted benchmark functions;
    shape holds s (step), w (window), S (fixed qubits), E (end qubits), h = S + E, b (bit-dependent extra width),
    i0, K0, the main term and D0 of the bound, and the bound"""
    funcs = {}
    for line in open(os.path.join(ROOT, "data", "benchmark", "pairs.jsonl")):
        p = json.loads(line)
        for side in ("ref", "sol"):
            funcs[(p[side]["file"], p[side]["function"])] = 1
    seen = set()
    for file, fn in sorted(funcs):
        d = analyze(os.path.join(ROOT, file), fn)
        if d["status"] != "ladder":
            continue
        fam = family(d)
        e = fam.extra
        for sg in fam.segments:
            if not isinstance(sg, Loop):
                continue
            key = repr((sg, e))
            if key in seen:
                continue
            seen.add(key)
            s0 = mirror(sg, e) if sg.step < 0 else sg
            refs = [q for g in s0.body for q in g.qubits] + [r for g in s0.body for r, _ in g.cond]
            offs = [k for kind, k in refs if kind == "i"]
            lo, w = min(offs), max(offs) - min(offs) + 1
            i0, dd = s0.start[1] + lo, s0.stop[1] + lo - e
            K = 1 - dd - w
            st = s0.step
            S = len({k for kind, k in refs if kind == "c"})
            E = len({k for kind, k in refs if kind == "n"})
            b = (w - 1 + S) if any(g.cond for g in s0.body) else 0
            h = S + E
            main = 4 ** h * 2 ** b * ((4 ** w - 1) // 3 + max(0, st - w))
            D0 = i0 * 8 ** S * 2 ** E + 4 ** E * (st * (K + 1) + st * (st - 1) // 2)
            yield fn, file, fam, sg, dict(s=st, w=w, S=S, E=E, h=h, b=b, i0=i0, K0=K, main_term=main, D0=D0,
                                          bound=main + D0)
