from qiskit import QuantumCircuit

def or_oracle(qc: QuantumCircuit, x: list, y: int) -> QuantumCircuit:
    """
    OR oracle that flips target if any query qubit is |1>.

    Args:
        qc: Quantum circuit to modify
        x: List of query qubit indices
        y: Target qubit index

    Returns:
        Modified quantum circuit
    """
    # OR can be implemented using De Morgan's law: A ∨ B = ¬(¬A ∧ ¬B)
    # Flip all query qubits
    for qubit in x:
        qc.x(qubit)

    # Apply multi-controlled X (AND of all flipped inputs)
    qc.mcx(x, y)

    # Flip target (to get NOT of the AND)
    qc.x(y)

    # Flip query qubits back
    for qubit in x:
        qc.x(qubit)


def solve(n):
    qc = QuantumCircuit(n + 1)
    or_oracle(qc, list(range(0, n)), n)
    return qc
