from qiskit import QuantumCircuit
import numpy as np

def binary_fraction_quantum_in_place(qc: QuantumCircuit, register: list) -> QuantumCircuit:
    """
    Apply in-place binary fraction transformation to the first qubit controlled by the rest.

    Args:
        qc: Quantum circuit to modify
        register: List of qubit indices in the register

    Returns:
        Modified quantum circuit
    """
    if len(register) == 0:
        return qc

    # Apply Hadamard to first qubit (equivalent to 1-qubit QFT)
    qc.h(register[0])

    # Apply controlled rotations from remaining qubits
    for ind in range(1, len(register)):
        # Apply controlled phase 2π / 2^(ind+1)
        phase = 2.0 * np.pi / (2**(ind + 1))
        qc.cp(phase, register[ind], register[0])

    return qc


def solve(n):
    qc = QuantumCircuit(n)
    binary_fraction_quantum_in_place(qc, list(range(0, n)))
    return qc
