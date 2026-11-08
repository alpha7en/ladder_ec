from qiskit import QuantumCircuit

def arbitrary_adder_in_place(qc: QuantumCircuit, a: list, b: list, carry: int) -> QuantumCircuit:
    n = len(a)
    max_qubit = max(max(a), max(b), carry)
    ancillas = [max_qubit + 1 + i for i in range(n)]

    # Forward: compute carries
    qc.ccx(a[0], b[0], ancillas[0])
    for i in range(1, n):
        qc.ccx(a[i], b[i], ancillas[i])
        qc.ccx(a[i], ancillas[i-1], ancillas[i])
        qc.ccx(b[i], ancillas[i-1], ancillas[i])
    qc.cx(ancillas[n-1], carry)

    # Backward: compute sums and clean
    for i in range(n-1, 0, -1):
        qc.ccx(b[i], ancillas[i-1], ancillas[i])
        qc.ccx(a[i], ancillas[i-1], ancillas[i])
        qc.ccx(a[i], b[i], ancillas[i])
        qc.cx(ancillas[i-1], b[i])
        qc.cx(a[i], b[i])

    qc.ccx(a[0], b[0], ancillas[0])
    qc.cx(a[0], b[0])

    return qc


def solve(n):
    qc = QuantumCircuit(2 * n + 1)
    arbitrary_adder_in_place(qc, list(range(0, n)), list(range(n, 2 * n)), 2 * n)
    return qc
