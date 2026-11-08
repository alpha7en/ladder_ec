"""Qiskit generators for the comparison with AutoQ-Para (Table 9): the same circuits as the GHZ models of
autoq_models.cpp and as Table 1 of the AutoQ-Para paper (Abdulla et al., 2025: BV, adder, syndrome, Grover,
Heisenberg Hamiltonian). The parameter n is the number of qubits of the circuit. The circuits of that table are
written from the descriptions in its Sec. 4.1-4.5, not translated from the AutoQ-Para automata."""
from qiskit import QuantumCircuit
import math


# ---------------- GHZ (as in autoq_models.cpp)
def ghz_ladder(n):
    c = QuantumCircuit(n); c.h(0)
    for i in range(n - 1):
        c.cx(i, i + 1)
    return c


def ghz_elevator(n):
    c = QuantumCircuit(n); c.h(0)
    for i in range(1, n):
        c.cx(0, i)
    return c


def ghz_ladder_noh(n):
    c = QuantumCircuit(n)
    for i in range(n - 1):
        c.cx(i, i + 1)
    return c


def ghz_ladder_rev(n):
    c = QuantumCircuit(n); c.h(0)
    for i in range(n - 1):
        c.cx(i + 1, i)
    return c


def cx_ladder(n):
    c = QuantumCircuit(n)
    for i in range(n - 1):
        c.cx(i, i + 1)
    return c


def cx_ladder_czh(n):
    c = QuantumCircuit(n)
    for i in range(n - 1):
        c.h(i + 1); c.cz(i, i + 1); c.h(i + 1)
    return c


def cx_ladder_rev(n):
    c = QuantumCircuit(n)
    for i in range(n - 1):
        c.cx(i + 1, i)
    return c


def cx_ladder_rev_hcxh(n):
    """model check of autoq_models.cpp: reversed CX as (H x H) CX (H x H)"""
    c = QuantumCircuit(n)
    for i in range(n - 1):
        c.h(i); c.h(i + 1); c.cx(i, i + 1); c.h(i); c.h(i + 1)
    return c


def cx_elevator(n):
    c = QuantumCircuit(n)
    for i in range(1, n):
        c.cx(0, i)
    return c


# ---------------- BV with hidden string (10)^* (Sec. 4.1 of AutoQ-Para): n - 1 data qubits, ancilla last
def bv10(n):
    k = n - 1
    c = QuantumCircuit(n); c.x(k)
    for i in range(n):
        c.h(i)
    for i in range(0, k, 2):
        c.cx(i, k)
    for i in range(k):
        c.h(i)
    return c


def bv10_post(n):
    """postcondition as a circuit: |s>|-> with s = 1010..."""
    k = n - 1
    c = QuantumCircuit(n)
    for i in range(0, k, 2):
        c.x(i)
    c.x(k); c.h(k)
    return c


def bv10_mut(n):
    """mutant: the oracle misses the first 1 bit of the hidden string"""
    k = n - 1
    c = QuantumCircuit(n); c.x(k)
    for i in range(n):
        c.h(i)
    for i in range(2, k, 2):
        c.cx(i, k)
    for i in range(k):
        c.h(i)
    return c


def bv_bits(n, a):
    """BV with any hidden string a (the bits of the task), to compare generality"""
    k = n - 1
    c = QuantumCircuit(n); c.x(k)
    for i in range(n):
        c.h(i)
    for i in range(k):
        if (a >> i) & 1:
            c.cx(i, k)
    for i in range(k):
        c.h(i)
    return c


def bv_bits_post(n, a):
    k = n - 1
    c = QuantumCircuit(n)
    for i in range(k):
        if (a >> i) & 1:
            c.x(i)
    c.x(k); c.h(k)
    return c


# ---------------- Cuccaro adder (Sec. 4.2 of AutoQ-Para, Fig. 10): c0, b0, a0, b1, a1, ..., ancilla z last
def _maj(c, x, y, z):
    c.cx(z, y); c.cx(z, x); c.ccx(x, y, z)


def _uma2(c, x, y, z):
    c.ccx(x, y, z); c.cx(z, x); c.cx(x, y)


def _uma3(c, x, y, z):
    c.x(y); c.cx(x, y); c.ccx(x, y, z); c.x(y); c.cx(z, x); c.cx(z, y)


def adder_uma2(n):
    c = QuantumCircuit(n)
    for i in range(0, n - 3, 2):
        _maj(c, i, i + 1, i + 2)
    c.cx(n - 2, n - 1)
    for i in range(n - 4, -1, -2):
        _uma2(c, i, i + 1, i + 2)
    return c


def adder_uma3(n):
    c = QuantumCircuit(n)
    for i in range(0, n - 3, 2):
        _maj(c, i, i + 1, i + 2)
    c.cx(n - 2, n - 1)
    for i in range(n - 4, -1, -2):
        _uma3(c, i, i + 1, i + 2)
    return c


def adder_mut(n):
    """mutant: the last UMA loses CX(z, x)"""
    c = QuantumCircuit(n)
    for i in range(0, n - 3, 2):
        _maj(c, i, i + 1, i + 2)
    c.cx(n - 2, n - 1)
    for i in range(n - 4, 1, -2):
        _uma2(c, i, i + 1, i + 2)
    c.ccx(0, 1, 2); c.cx(0, 1)
    return c


# ---------------- syndrome extraction of the repetition code (Sec. 4.3, Fig. 11): x_1, a_1, x_2, a_2, ...
def syndrome(n):
    """a_i ^= x_{i-1} ^ x_i, i >= 2; qubits x_i = 2(i-1), a_i = 2(i-1)+1"""
    c = QuantumCircuit(n)
    for j in range(3, n, 2):
        c.cx(j - 3, j); c.cx(j - 1, j)
    return c


def syndrome_two_passes(n):
    """the same CX, first all CX(x_i, a_i), then all CX(x_{i-1}, a_i)"""
    c = QuantumCircuit(n)
    for j in range(3, n, 2):
        c.cx(j - 1, j)
    for j in range(3, n, 2):
        c.cx(j - 3, j)
    return c


def syndrome_mut(n):
    """mutant: the last a_i has no CX from x_{i-1}"""
    c = QuantumCircuit(n)
    for j in range(3, n - 2, 2):
        c.cx(j - 3, j); c.cx(j - 1, j)
    c.cx(n - 2, n - 1)
    return c


# ---------------- Grover, one iteration (Sec. 4.4, Fig. 12a): multi-controlled gates
def grover_ideal(n):
    c = QuantumCircuit(n)
    for i in range(n):
        c.h(i)
    c.mcx(list(range(n - 1)), n - 1)
    for i in range(n):
        c.h(i)
    return c


# ---------------- Heisenberg Hamiltonian, one Trotter step (Sec. 4.5, Figs. 13e and 14), 2 delta = pi/4
def heis_naive(n):
    """for every pair: XX = (H x H) Rzz (H x H), YY = (S x S)(H x H) Rzz (H x H)(Sdg x Sdg), ZZ = Rzz"""
    c = QuantumCircuit(n)
    for j in range(n - 1):
        c.h(j); c.h(j + 1); c.rzz(math.pi / 4, j, j + 1); c.h(j); c.h(j + 1)
    for j in range(n - 1):
        c.sdg(j); c.sdg(j + 1); c.h(j); c.h(j + 1); c.rzz(math.pi / 4, j, j + 1); c.h(j); c.h(j + 1); c.s(j); c.s(j + 1)
    for j in range(n - 1):
        c.rzz(math.pi / 4, j, j + 1)
    return c


def heis_opt(n):
    """after cancelling H H = I and S Sdg = I between neighbouring boxes: layers of single-qubit gates and three Rzz ladders"""
    c = QuantumCircuit(n)
    for j in range(n):
        c.h(j)
    for j in range(n - 1):
        c.rzz(math.pi / 4, j, j + 1)
    for j in range(n):
        c.h(j); c.sdg(j); c.h(j)
    for j in range(n - 1):
        c.rzz(math.pi / 4, j, j + 1)
    for j in range(n):
        c.h(j); c.s(j)
    for j in range(n - 1):
        c.rzz(math.pi / 4, j, j + 1)
    return c


def heis_mut(n):
    """mutant: the second layer loses S on the last qubit"""
    c = QuantumCircuit(n)
    for j in range(n):
        c.h(j)
    for j in range(n - 1):
        c.rzz(math.pi / 4, j, j + 1)
    for j in range(n):
        c.h(j); c.sdg(j); c.h(j)
    for j in range(n - 1):
        c.rzz(math.pi / 4, j, j + 1)
    for j in range(n - 1):
        c.h(j); c.s(j)
    c.h(n - 1)
    for j in range(n - 1):
        c.rzz(math.pi / 4, j, j + 1)
    return c
