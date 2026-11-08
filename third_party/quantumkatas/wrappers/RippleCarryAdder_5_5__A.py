from qiskit import QuantumCircuit

def in_place_subtractor_modulo_n(qc: QuantumCircuit, a: list, b: list) -> QuantumCircuit:
    # Subtraction is the adjoint of addition in modular arithmetic
    n = len(a)
    max_qubit = max(max(a), max(b))
    ancilla = max_qubit + 1
    carries = [ancilla] + a

    # Reverse UMA phase (forward direction of adjoint)
    for i in range(n):
        qc.cx(carries[i], b[i])
        qc.cx(a[i], carries[i])
        qc.ccx(b[i], carries[i], a[i])

    # Reverse Majority phase
    for i in range(n-1, -1, -1):
        qc.ccx(b[i], carries[i], a[i])
        qc.cx(a[i], carries[i])
        qc.cx(a[i], b[i])

    return qc


def solve(n):
    qc = QuantumCircuit(2 * n)
    in_place_subtractor_modulo_n(qc, list(range(0, n)), list(range(n, 2 * n)))
    return qc
