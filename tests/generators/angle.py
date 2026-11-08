import math
from qiskit import QuantumCircuit


def rz_residue(n):
    qc = QuantumCircuit(n)
    for i in range(n):
        qc.rz(math.pi / 4 + 0.1, i)
    return qc


def rz_exact(n):
    qc = QuantumCircuit(n)
    for i in range(n):
        qc.rz(math.pi / 4, i)
    return qc
