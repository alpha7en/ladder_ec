from qiskit import QuantumCircuit
from typing import List
import numpy as np

def quantum_fourier_pattern(qc: QuantumCircuit, qubits: List[int]) -> QuantumCircuit:
    """
    Apply the Quantum Fourier Transform to the given qubits (little-endian).

    Args:
        qc: Quantum circuit to modify
        qubits: List of qubit indices for QFT (qubits[0] is the least-significant input bit)

    Returns:
        Modified quantum circuit
    """
    n = len(qubits)

    # Standard QFT: process high-order bit first, apply H then controlled phases
    # from subsequent (lower-order) bits.
    for i in reversed(range(n)):
        qc.h(qubits[i])
        for j in range(i):
            angle = np.pi / (2 ** (i - j))
            qc.cp(angle, qubits[j], qubits[i])

    # Final swaps to match the conventional QFT output ordering.
    for k in range(n // 2):
        qc.swap(qubits[k], qubits[n - 1 - k])

    return qc


def solve(n):
    qc = QuantumCircuit(n)
    quantum_fourier_pattern(qc, list(range(0, n)))
    return qc
