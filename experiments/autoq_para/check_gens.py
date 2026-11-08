"""Ground truth of the AutoQ-Para pairs by direct computation: Qiskit matrices at n = 3..9 (for tasks with bits,
all a). Usage: python -I experiments/autoq_para/check_gens.py -> data/runs/autoq_para/truth.json"""
import itertools, json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import numpy as np                                        # noqa: E402
from qiskit.quantum_info import Operator, Statevector    # noqa: E402
import gens                                               # noqa: E402
from pairs import PAIRS                                   # noqa: E402

out = {}
for name, r, s, rel, exp, *_ in PAIRS:
    fr, fs = getattr(gens, r), getattr(gens, s)
    bits = fr.__code__.co_argcount == 2
    ok, bad = True, []
    for n in range(3, 10):
        for a in (itertools.product((0, 1), repeat=n) if bits else [None]):
            args = (n, sum(b << i for i, b in enumerate(a))) if bits else (n,)
            A, B = fr(*args), fs(*args)
            if rel == "E1":
                eq = np.allclose(Operator(A).data, Operator(B).data, atol=1e-9)
            else:
                eq = np.allclose(Statevector(A).data, Statevector(B).data, atol=1e-9)
            if not eq:
                ok = False
                bad.append(n)
                break
    out[name] = dict(truth=ok, expected=exp, first_bad_n=bad[:1])
    print(name, ok, "expected", exp, bad[:1])
json.dump(out, open(os.path.join(HERE, "..", "..", "data", "runs", "autoq_para", "truth.json"), "w"), indent=1,
          ensure_ascii=False)
