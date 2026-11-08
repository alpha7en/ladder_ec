from qiskit import QuantumCircuit
import numpy as np
from typing import List

def wstate_arbitrary(qc: QuantumCircuit, qs: List[int]) -> QuantumCircuit:
    n = len(qs)

    if n == 1:
        qc.x(qs[0])
        return qc

    qc.x(qs[0])
    for i in range(n - 1):
        theta = 2 * np.arccos(np.sqrt(1 / (n - i)))
        qc.cry(theta, qs[i], qs[i + 1])
        qc.cx(qs[i + 1], qs[i])


def solve(n):
    qc = QuantumCircuit(n)
    wstate_arbitrary(qc, list(range(0, n)))
    return qc
