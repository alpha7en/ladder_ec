from qiskit import QuantumCircuit

def hadamard_transform(qc: QuantumCircuit, register: list[int]) -> QuantumCircuit:
    """
    Apply Hadamard gate to all qubits in the register.

    Args:
        qc: Quantum circuit to modify
        register: List of qubit indices

    Returns:
        Modified quantum circuit
    """
    # Apply Hadamard gate to each qubit in the register
    for qubit in register:
        qc.h(qubit)

    return qc


def solve(n):
    qc = QuantumCircuit(n)
    hadamard_transform(qc, list(range(0, n)))
    return qc
