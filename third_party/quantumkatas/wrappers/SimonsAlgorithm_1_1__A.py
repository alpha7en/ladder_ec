from qiskit import QuantumCircuit

def oracle_count_bits(qc: QuantumCircuit, x_register: list[int], y: int) -> QuantumCircuit:
    """
    Oracle that computes the parity of all input bits.

    Args:
        qc: Quantum circuit to modify
        x_register: List of input qubit indices
        y: Output qubit index

    Returns:
        Modified quantum circuit
    """
    # Apply CNOT from each input qubit to output qubit
    # This computes y ⊕ (x₀ ⊕ x₁ ⊕ ... ⊕ xₙ₋₁)
    for x_qubit in x_register:
        qc.cx(x_qubit, y)
    return qc


def solve(n):
    qc = QuantumCircuit(n + 1)
    oracle_count_bits(qc, list(range(0, n)), n)
    return qc
