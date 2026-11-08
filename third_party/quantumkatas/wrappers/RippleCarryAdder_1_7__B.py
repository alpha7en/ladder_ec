from qiskit import QuantumCircuit

def arbitrary_adder(qc: QuantumCircuit, a: list, b: list, sum_qubits: list, carry: int) -> QuantumCircuit:
    n = len(a)
    if n == 1:
        qc.cx(a[0], sum_qubits[0])
        qc.cx(b[0], sum_qubits[0])
        qc.ccx(a[0], b[0], carry)
        return qc

    # Need n-1 ancilla qubits for internal carries
    max_qubit = max(max(a), max(b), max(sum_qubits), carry)
    ancillas = [max_qubit + 1 + i for i in range(n-1)]

    # Forward: compute carries
    qc.ccx(a[0], b[0], ancillas[0])
    for i in range(1, n-1):
        qc.ccx(a[i], b[i], ancillas[i])
        qc.ccx(a[i], ancillas[i-1], ancillas[i])
        qc.ccx(b[i], ancillas[i-1], ancillas[i])

    # Final carry
    qc.ccx(a[n-1], b[n-1], carry)
    qc.ccx(a[n-1], ancillas[n-2], carry)
    qc.ccx(b[n-1], ancillas[n-2], carry)

    # Compute sums
    qc.cx(a[0], sum_qubits[0])
    qc.cx(b[0], sum_qubits[0])
    for i in range(1, n):
        qc.cx(a[i], sum_qubits[i])
        qc.cx(b[i], sum_qubits[i])
        qc.cx(ancillas[i-1], sum_qubits[i])

    # Cleanup carries (reverse order)
    for i in range(n-2, 0, -1):
        qc.ccx(b[i], ancillas[i-1], ancillas[i])
        qc.ccx(a[i], ancillas[i-1], ancillas[i])
        qc.ccx(a[i], b[i], ancillas[i])
    qc.ccx(a[0], b[0], ancillas[0])


def solve(n):
    qc = QuantumCircuit(3 * n + 1)
    arbitrary_adder(qc, list(range(0, n)), list(range(n, 2 * n)), list(range(2 * n, 3 * n)), 3 * n)
    return qc
