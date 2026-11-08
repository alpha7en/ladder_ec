from qiskit import QuantumCircuit

def oracle_product_function(qc: QuantumCircuit, x: list[int], y: int, r: list[int]) -> QuantumCircuit:
    """
    Oracle for scalar product function f(x) = Σᵢ rᵢ xᵢ mod 2.

    Args:
        qc: Quantum circuit to modify
        x: List of input qubit indices
        y: Output qubit index
        r: Bit vector defining the function

    Returns:
        Modified quantum circuit
    """
    assert len(x) == len(r), "Arrays should have the same length"

    # For each bit in r that is 1, apply CNOT with corresponding input qubit as control
    for i, bit in enumerate(r):
        if bit == 1:
            qc.cx(x[i], y)

    return qc


def solve(n, a):
    qc = QuantumCircuit(n + 1)
    oracle_product_function(qc, list(range(0, n)), n, a)
    return qc
