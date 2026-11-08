from qiskit import QuantumCircuit

def oracle_bitwise_right_shift(qc: QuantumCircuit, x_register: list[int], y_register: list[int]) -> QuantumCircuit:
    """
    Oracle that applies bitwise right shift function.

    Args:
        qc: Quantum circuit to modify
        x_register: List of input qubit indices
        y_register: List of output qubit indices (same length as x_register)

    Returns:
        Modified quantum circuit
    """
    n = len(x_register)
    # Right shift: y[i] ⊕= x[i-1] for i > 0
    # y[0] remains unchanged (no bit to shift from)
    for i in range(1, n):
        qc.cx(x_register[i-1], y_register[i])
    return qc


def solve(n):
    qc = QuantumCircuit(2 * n)
    oracle_bitwise_right_shift(qc, list(range(0, n)), list(range(n, 2 * n)))
    return qc
