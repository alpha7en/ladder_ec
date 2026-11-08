from qiskit import QuantumCircuit
from typing import List
import numpy as np

def inverse_qft(qc: QuantumCircuit, register: List[int]) -> QuantumCircuit:
    """
    Apply inverse quantum Fourier transform.

    Args:
        qc: Quantum circuit to modify
        register: List of qubit indices

    Returns:
        Modified quantum circuit
    """
    n = len(register)

    # Inverse QFT is the adjoint of QFT
    # First reverse the qubit order (reverse the bit reversal step)
    for i in range(n // 2):
        qc.swap(register[i], register[n - 1 - i])

    # Apply inverse rotations in reverse order
    for i in range(n):
        # Apply controlled rotations in reverse
        for j in range(i):
            # Phase rotation R_k† = R_{-k}
            k = i - j + 1
            qc.cp(-2 * np.pi / (2**k), register[j], register[i])

        # Inverse Hadamard at the end
        qc.h(register[i])

    return qc


def solve(n):
    qc = QuantumCircuit(n)
    inverse_qft(qc, list(range(0, n)))
    return qc
