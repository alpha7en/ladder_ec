"""Writes gens.py of the step-4 experiment (Table 8): every function in full (the translator does not accept a call
of a helper that builds the circuit). Usage: python -I experiments/step4/make_gens.py"""
import os

REF = {
 "cx": ("(n)", "    qc = QuantumCircuit(n)\n    for i in range(0, n - 1):\n        qc.cx(i, i + 1)\n",
        "        qc.cx(i, i + 1)\n        qc.cx(i, i + 1)\n"),
 "rzz": ("(n)", "    qc = QuantumCircuit(n)\n    for i in range(n):\n        qc.h(i)\n    for i in range(0, n - 1):\n        qc.rzz(math.pi / 4, i, i + 1)\n",
         "        qc.rzz(math.pi / 8, i, i + 1)\n        qc.rzz(-math.pi / 8, i, i + 1)\n"),
 "bv": ("(n, a)", "    qc = QuantumCircuit(n)\n    qc.x(n - 1)\n    for i in range(n):\n        qc.h(i)\n    for i in range(n - 1):\n        if (a >> i) & 1:\n            qc.cx(i, n - 1)\n    for i in range(n - 1):\n        qc.h(i)\n",
        "        if (a >> i) & 1:\n            qc.cx(i, n - 1)\n            qc.cx(i, n - 1)\n"),
}
KS = {"cx": (5, 8, 11, 14, 20), "rzz": (5, 8, 11, 14, 20), "bv": (5, 8, 11)}
out = ['"""Step-4 experiment (Table 8): correct solutions with a constant K. Solution = reference + an extra loop from K\n'
       'whose body cancels (a block and its inverse), like the mutation C6_K of the benchmark, but without an error.\n'
       'Written by make_gens.py; do not edit by hand."""', "import math", "from qiskit import QuantumCircuit", ""]
for fam, (sig, body, extra) in REF.items():
    out += ["", f"def {fam}_ref{sig}:", body + "    return qc", ""]
    for K in KS[fam]:
        out += ["", f"def {fam}_k{K}{sig}:", body + f"    for i in range({K}, n - 1):\n" + extra + "    return qc", ""]
open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "gens.py"), "w").write("\n".join(out) + "\n")
