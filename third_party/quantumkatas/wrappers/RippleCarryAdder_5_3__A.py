from qiskit import QuantumCircuit
import numpy as np

def subtractor_modulo_n(qc: QuantumCircuit, a: list, b: list, diff: list) -> QuantumCircuit:
    n = len(a)
    # Copy b to diff
    for i in range(n):
        qc.cx(b[i], diff[i])
    # QFT on diff (no swaps, MSB to LSB)
    for j in range(n-1, -1, -1):
        qc.h(diff[j])
        for k in range(j-1, -1, -1):
            qc.cp(np.pi / 2**(j-k), diff[k], diff[j])
    # Subtract a in Fourier space (negative phase rotations)
    for j in range(n):
        for i in range(j+1):
            qc.cp(-2 * np.pi / 2**(j - i + 1), a[i], diff[j])
    # Inverse QFT on diff
    for j in range(n):
        for k in range(j-1, -1, -1):
            qc.cp(-np.pi / 2**(j-k), diff[k], diff[j])
        qc.h(diff[j])
    return qc


def solve(n):
    qc = QuantumCircuit(3 * n)
    subtractor_modulo_n(qc, list(range(0, n)), list(range(n, 2 * n)), list(range(2 * n, 3 * n)))
    return qc
