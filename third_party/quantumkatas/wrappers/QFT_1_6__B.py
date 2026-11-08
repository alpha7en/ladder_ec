from qiskit import QuantumCircuit

def reverse_register(qc: QuantumCircuit, register: list) -> QuantumCircuit:
    """
    Reverse the order of qubits in the register.

    Args:
        qc: Quantum circuit to modify
        register: List of qubit indices to reverse

    Returns:
        Modified quantum circuit
    """
    n = len(register)
    # Swap qubits symmetrically from the ends
    for i in range(n // 2):
        qc.swap(register[i], register[n - 1 - i])


def solve(n):
    qc = QuantumCircuit(n)
    reverse_register(qc, list(range(0, n)))
    return qc
