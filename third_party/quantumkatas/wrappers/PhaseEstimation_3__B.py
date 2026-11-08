from qiskit import QuantumCircuit
from typing import List
import numpy as np

def inverse_qft(qc: QuantumCircuit, qubits: List[int]) -> QuantumCircuit:
    """
    Apply inverse Quantum Fourier Transform.

    Args:
        qc: Quantum circuit to modify
        qubits: List of qubit indices for inverse QFT

    Returns:
        Modified quantum circuit with inverse QFT applied
    """
    n = len(qubits)

    # Swap qubits first (reverse of QFT)
    for i in range(n // 2):
        qc.swap(qubits[i], qubits[n - 1 - i])

    # Apply inverse rotations in reverse order
    for i in range(n - 1, -1, -1):
        # Inverse controlled rotations
        for j in range(n - 1, i, -1):
            angle = -np.pi / (2 ** (j - i))
            qc.cp(angle, qubits[j], qubits[i])

        # Hadamard on qubit i
        qc.h(qubits[i])


def solve(n):
    qc = QuantumCircuit(n)
    inverse_qft(qc, list(range(0, n)))
    return qc
