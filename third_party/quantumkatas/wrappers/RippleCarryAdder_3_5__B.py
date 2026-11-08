from qiskit import QuantumCircuit

def arbitrary_maj_uma_adder(qc: QuantumCircuit, a: list, b: list, carry: int) -> QuantumCircuit:
    n = len(a)
    ancilla = max(max(a), max(b), carry) + 1
    carries = [ancilla] + a

    # Majority phase
    for i in range(n):
        qc.cx(a[i], b[i])
        qc.cx(a[i], carries[i])
        qc.ccx(b[i], carries[i], a[i])

    # Save final carry
    qc.cx(carries[n], carry)

    # UMA phase (reverse order)
    for i in range(n-1, -1, -1):
        qc.ccx(b[i], carries[i], a[i])
        qc.cx(a[i], carries[i])
        qc.cx(carries[i], b[i])


def solve(n):
    qc = QuantumCircuit(2 * n + 1)
    arbitrary_maj_uma_adder(qc, list(range(0, n)), list(range(n, 2 * n)), 2 * n)
    return qc
