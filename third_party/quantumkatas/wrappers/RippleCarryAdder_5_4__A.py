from qiskit import QuantumCircuit

def in_place_adder_modulo_n(qc: QuantumCircuit, a: list, b: list) -> QuantumCircuit:
    n = len(a)
    max_qubit = max(max(a), max(b))
    ancilla = max_qubit + 1
    carries = [ancilla] + a

    # Majority phase (don't compute final carry for modular arithmetic)
    for i in range(n):
        qc.cx(a[i], b[i])
        qc.cx(a[i], carries[i])
        qc.ccx(b[i], carries[i], a[i])

    # UMA phase (reverse order)
    for i in range(n-1, -1, -1):
        qc.ccx(b[i], carries[i], a[i])
        qc.cx(a[i], carries[i])
        qc.cx(carries[i], b[i])

    return qc


def solve(n):
    qc = QuantumCircuit(2 * n)
    in_place_adder_modulo_n(qc, list(range(0, n)), list(range(n, 2 * n)))
    return qc
