from qiskit import QuantumCircuit

def prepare_alternating_state(qc: QuantumCircuit, register: list) -> QuantumCircuit:
    """
    Prepare alternating amplitude state.
    We just need to create state 1/sqrt(N) * sum (-1)^k |k>.
    which is |-> ... |->.
    Wait, |->|-> is 1 -1 -1 1 (+ - - +).
    (+ - + -) is |+> ... |+> |->. (LSB is -).
    So H all, Z(register[0]).
    """
    n = len(register)

    for i in range(n):
        qc.h(register[i])

    qc.z(register[0])

    return qc


def solve(n):
    qc = QuantumCircuit(n)
    prepare_alternating_state(qc, list(range(0, n)))
    return qc
