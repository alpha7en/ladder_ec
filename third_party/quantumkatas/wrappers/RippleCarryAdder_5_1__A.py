from qiskit import QuantumCircuit

def adder_modulo_n(qc: QuantumCircuit, a: list, b: list, sum_qubits: list) -> QuantumCircuit:
    n = len(a)

    # Use sum qubits to store intermediate carries
    qc.ccx(a[0], b[0], sum_qubits[0])
    for i in range(1, n):
        qc.ccx(a[i], b[i], sum_qubits[i])
        qc.ccx(a[i], sum_qubits[i-1], sum_qubits[i])
        qc.ccx(b[i], sum_qubits[i-1], sum_qubits[i])

    # Compute sums and clean up carries
    for i in range(n-1, 0, -1):
        qc.ccx(b[i], sum_qubits[i-1], sum_qubits[i])
        qc.ccx(a[i], sum_qubits[i-1], sum_qubits[i])
        qc.ccx(a[i], b[i], sum_qubits[i])
        qc.cx(a[i], sum_qubits[i])
        qc.cx(b[i], sum_qubits[i])
        qc.cx(sum_qubits[i-1], sum_qubits[i])

    qc.ccx(a[0], b[0], sum_qubits[0])
    qc.cx(a[0], sum_qubits[0])
    qc.cx(b[0], sum_qubits[0])

    return qc


def solve(n):
    qc = QuantumCircuit(3 * n)
    adder_modulo_n(qc, list(range(0, n)), list(range(n, 2 * n)), list(range(2 * n, 3 * n)))
    return qc
