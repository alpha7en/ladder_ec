from qiskit import QuantumCircuit
from typing import List
import numpy as np

def qft(qc: QuantumCircuit, qubits: List[int]) -> QuantumCircuit:
    """
    Apply Quantum Fourier Transform.

    Args:
        qc: Quantum circuit to modify
        qubits: List of qubit indices for QFT

    Returns:
        Modified quantum circuit with QFT applied
    """
    n = len(qubits)

    for i in range(n):
        # Hadamard on qubit i
        qc.h(qubits[i])

        # Controlled rotations
        for j in range(i + 1, n):
            angle = np.pi / (2 ** (j - i))
            qc.cp(angle, qubits[j], qubits[i])

    # Swap qubits for correct ordering
    for i in range(n // 2):
        qc.swap(qubits[i], qubits[n - 1 - i])


def solve(n):
    qc = QuantumCircuit(n)
    qft(qc, list(range(0, n)))
    return qc
