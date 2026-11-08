from qiskit import QuantumCircuit
import numpy as np

def binary_fraction_quantum_in_place_helper(qc: QuantumCircuit, register: list) -> QuantumCircuit:
    """Helper function for in-place binary fraction transformation."""
    if len(register) == 0:
        return qc

    # Apply Hadamard to first qubit
    qc.h(register[0])

    # Apply controlled rotations from remaining qubits
    for ind in range(1, len(register)):
        phase = 2.0 * np.pi / (2**(ind + 1))
        qc.cp(phase, register[ind], register[0])

    return qc

def reverse_register_helper(qc: QuantumCircuit, register: list) -> QuantumCircuit:
    """Helper function to reverse register."""
    n = len(register)
    for i in range(n // 2):
        qc.swap(register[i], register[n - 1 - i])
    return qc

def quantum_fourier_transform(qc: QuantumCircuit, register: list) -> QuantumCircuit:
    """
    Apply Quantum Fourier Transform to the register.

    Args:
        qc: Quantum circuit to modify
        register: List of qubit indices to transform

    Returns:
        Modified quantum circuit
    """
    n = len(register)

    # Apply binary fraction transformation to each qubit in sequence
    for i in range(n):
        qc = binary_fraction_quantum_in_place_helper(qc, register[i:])

    # Reverse the order of qubits
    qc = reverse_register_helper(qc, register)

    return qc


def solve(n):
    qc = QuantumCircuit(n)
    quantum_fourier_transform(qc, list(range(0, n)))
    return qc
