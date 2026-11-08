from qiskit import QuantumCircuit

def oracle_operator_output(qc: QuantumCircuit, x: list[int], y: int, A: list[int]) -> QuantumCircuit:
    """
    Oracle for linear operator function f(x) = A·x.

    Args:
        qc: Quantum circuit to modify
        x: List of input qubit indices
        y: Output qubit index
        A: Binary matrix as list (1xN)

    Returns:
        Modified quantum circuit
    """
    assert len(x) == len(A), "Arrays x and A should have the same length"

    # Linear operator A(x) = Σᵢ A[i] * x[i]
    for i, coeff in enumerate(A):
        if coeff == 1:
            qc.cx(x[i], y)

    return qc


def solve(n, a):
    qc = QuantumCircuit(n + 1)
    oracle_operator_output(qc, list(range(0, n)), n, a)
    return qc
