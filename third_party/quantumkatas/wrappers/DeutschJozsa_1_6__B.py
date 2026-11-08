from qiskit import QuantumCircuit

def oracle_product_with_negation_function(qc: QuantumCircuit, x: list[int], y: int, r: list[int]) -> QuantumCircuit:
    """
    Oracle for function f(x) = Σᵢ (rᵢ xᵢ + (1 - rᵢ)(1 - xᵢ)) mod 2.

    Args:
        qc: Quantum circuit to modify
        x: List of input qubit indices
        y: Output qubit index
        r: Bit vector defining the function

    Returns:
        Modified quantum circuit
    """
    assert len(x) == len(r), "Arrays should have the same length"

    # This function can be simplified:
    # f(x) = Σᵢ (rᵢ xᵢ + (1 - rᵢ)(1 - xᵢ))
    #      = Σᵢ (rᵢ xᵢ + (1 - rᵢ - xᵢ + rᵢ xᵢ))
    #      = Σᵢ (2 rᵢ xᵢ + 1 - rᵢ - xᵢ)
    #      = Σᵢ (1 - rᵢ - xᵢ) (since 2 rᵢ xᵢ ≡ 0 mod 2)
    #      = Σᵢ (1 - rᵢ) + Σᵢ (- xᵢ)
    #      = Σᵢ (1 - rᵢ) - Σᵢ xᵢ
    #      = |r| - Σᵢ rᵢ - Σᵢ xᵢ (where |r| is number of bits in r)
    #
    # But let's implement it more directly:
    # For each i: if r[i] == 1, contribute x[i]; if r[i] == 0, contribute (1 - x[i])

    # First, add contribution of (1 - r[i]) for each r[i] == 0
    for i, bit in enumerate(r):
        if bit == 0:
            # Contribute 1, which means flip the output
            qc.x(y)

    # Then, for each qubit, contribute based on r[i]
    for i, bit in enumerate(r):
        if bit == 1:
            # Contribute x[i] directly
            qc.cx(x[i], y)
        else:
            # Contribute (1 - x[i]) = 1 + x[i], but we already added 1 above, so subtract x[i]
            qc.cx(x[i], y)


def solve(n, a):
    qc = QuantumCircuit(n + 1)
    oracle_product_with_negation_function(qc, list(range(0, n)), n, a)
    return qc
