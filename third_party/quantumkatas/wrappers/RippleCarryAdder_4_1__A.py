from qiskit import QuantumCircuit

def subtractor(qc: QuantumCircuit, a: list, b: list, borrow_bit: int) -> QuantumCircuit:
    # Flip b to get one's complement
    for qubit in b:
        qc.x(qubit)

    # Add a + (~b) using existing adder (simulated with Maj-UMA)
    n = len(a)
    max_qubit = max(max(a), max(b), borrow_bit)
    ancilla = max_qubit + 1
    carries = [ancilla] + a

    # Majority phase
    for i in range(n):
        qc.cx(a[i], b[i])
        qc.cx(a[i], carries[i])
        qc.ccx(b[i], carries[i], a[i])

    qc.cx(carries[n], borrow_bit)

    # UMA phase
    for i in range(n-1, -1, -1):
        qc.ccx(b[i], carries[i], a[i])
        qc.cx(a[i], carries[i])
        qc.cx(carries[i], b[i])

    # Flip b back to get two's complement result
    for qubit in b:
        qc.x(qubit)

    return qc


def solve(n):
    qc = QuantumCircuit(2 * n + 1)
    subtractor(qc, list(range(0, n)), list(range(n, 2 * n)), 2 * n)
    return qc
