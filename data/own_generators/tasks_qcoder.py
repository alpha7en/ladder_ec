"""Parameterized tasks named after problems of the QCoder contests (the field source gives the problem id). The
problem statements are not included; each task has a one-line English title and a specification written by us.
All solutions are our own and were checked by simulation against the specification."""
from __future__ import annotations

import math
import random

import numpy as np
from qiskit import QuantumCircuit
from qiskit.circuit.library import HGate, RYGate, StatePreparation

from common import (Style, Task, bits_of, iqft_textbook, mcx, mcz_pattern,
                    perm_matrix, qft_textbook, state_spec, unitary_spec)


def _hn(n):
    H = np.array([[1, 1], [1, -1]]) / math.sqrt(2)
    M = np.array([[1]], dtype=complex)
    for _ in range(n):
        M = np.kron(H, M)
    return M


# =====================================================================
# Q01. (|0...0> - |1...1>)/sqrt(2) with depth at most 6
# =====================================================================

def ghzm_ladder(n):
    qc = QuantumCircuit(n)
    qc.h(0); qc.z(0)
    for i in range(n - 1):
        qc.cx(i, i + 1)
    return qc


def ghzm_tree(n):
    qc = QuantumCircuit(n)
    qc.x(0); qc.h(0)
    have = 1
    while have < n:
        step = min(have, n - have)
        for i in range(step):
            qc.cx(i, have + i)
        have += step
    return qc


def _ghzm_spec(n):
    v = np.zeros(2 ** n); v[0] = 1; v[-1] = -1
    return state_spec(v)


Q01 = Task(
    id="Q01_ghz_minus", title="(|0...0> - |1...1>)/sqrt(2) with depth at most 6",
    source="QPC002_A3 / A4 / A5", equiv="E5", ptype="P2", params=["n"], spec=_ghzm_spec,
    grid=lambda: [{"n": n} for n in range(2, 8)],
    styles=[
        Style("ladder", ghzm_ladder, "M2", [], ""),
        Style("tree", ghzm_tree, "M4-tree", ["log-depth", "index-2^t"],
              ""),
    ],
    depth_limit=lambda n: 6, qcoder_range=lambda: [{"n": n} for n in range(2, 16)],
)


# =====================================================================
# Q02. Uniform superposition H^n
# =====================================================================

def unif_loop(n):
    qc = QuantumCircuit(n)
    for i in range(n):
        qc.h(i)
    return qc


def unif_broadcast(n):
    qc = QuantumCircuit(n)
    qc.h(range(n))
    return qc


Q02 = Task(
    id="Q02_uniform", title="uniform superposition H^n", source="QPC001_A2 = QPC004_B1",
    equiv="E5", ptype="P2", params=["n"],
    spec=lambda n: state_spec(np.ones(2 ** n)), grid=lambda: [{"n": n} for n in range(1, 8)],
    styles=[
        Style("loop", unif_loop, "M2", [], "for i in range(n): h(i)"),
        Style("broadcast", unif_broadcast, "M2", ["broadcast"],
              ""),
    ],
)


# =====================================================================
# Q03. Parity oracle
# =====================================================================

def parity_star(n):
    qc = QuantumCircuit(n + 1)
    for i in range(n):
        qc.cx(i, n)
    return qc


def parity_ladder(n):
    qc = QuantumCircuit(n + 1)
    for i in range(n - 1):
        qc.cx(i, i + 1)
    qc.cx(n - 1, n)
    for i in reversed(range(n - 1)):
        qc.cx(i, i + 1)
    return qc


def parity_tree(n):
    """The parity is folded by a tree into qubit 0, copied to y, and the tree is unfolded."""
    qc = QuantumCircuit(n + 1)
    layers = []
    stride = 1
    while stride < n:
        layer = [(i + stride, i) for i in range(0, n, 2 * stride) if i + stride < n]
        layers.append(layer)
        stride *= 2
    for layer in layers:
        for c, t in layer:
            qc.cx(c, t)
    qc.cx(0, n)
    for layer in reversed(layers):
        for c, t in layer:
            qc.cx(c, t)
    return qc


def _parity_spec(n):
    return unitary_spec(perm_matrix(n + 1, lambda x: x ^ ((bin(x % 2 ** n).count("1") % 2) << n)))


Q03 = Task(
    id="Q03_parity", title="parity oracle", source="QPC001_B2",
    equiv="E2", ptype="P2", params=["n"], spec=_parity_spec,
    grid=lambda: [{"n": n} for n in range(1, 7)],
    styles=[
        Style("star", parity_star, "M3", [], ""),
        Style("ladder", parity_ladder, "M2", ["two-ladders"], ""),
        Style("tree", parity_tree, "M4-tree", ["log-depth", "index-2^t"], ""),
    ],
)


# =====================================================================
# Q04. Phase oracle x < L
# =====================================================================

def flipL_enumerate(n, L):
    """Enumerates all x < L, one multi-controlled Z each: of the order of L gates, up to 2^n."""
    qc = QuantumCircuit(n)
    for x in range(L):
        mcz_pattern(qc, range(n), bits_of(x, n))
    return qc


def flipL_compare(n, L):
    """x < L iff there is a bit t with L_t = 1, x_t = 0 and equal higher bits; these conditions are disjoint, so n phase gates suffice."""
    qc = QuantumCircuit(n)
    if L == 2 ** n:
        return qc
    for t in range(n):
        if (L >> t) & 1:
            qs = list(range(t, n))
            bits = [0] + [(L >> j) & 1 for j in range(t + 1, n)]
            mcz_pattern(qc, qs, bits)
    return qc


def _flipL_spec(n, L):
    return unitary_spec(np.diag([-1 if x < L else 1 for x in range(2 ** n)]).astype(complex))


def _flipL_grid():
    return [{"n": n, "L": L} for n in range(1, 5) for L in range(1, 2 ** n + 1)]


Q04 = Task(
    id="Q04_flip_below_L", title="phase oracle x < L", source="QPC001_B3 / B4",
    equiv="E2", ptype="P5", params=["n", "L"], spec=_flipL_spec, grid=_flipL_grid,
    styles=[
        Style("enumerate", flipL_enumerate, "nonM", ["if-data", "mcx", "count-L"],
              ""),
        Style("compare", flipL_compare, "M3*", ["if-data", "mcx"], ""),
    ],
    depth_limit=lambda n, L: 50,
    qcoder_range=lambda: [{"n": n, "L": L} for n in (8, 10) for L in (1, 2 ** n // 3, 2 ** n - 1)],
    note="",
)


# =====================================================================
# Q05. Uniform superposition of the first L basis states
# =====================================================================

def _ctrl_apply(qc, gate, target, ctrls):
    if not ctrls:
        qc.append(gate, [target])
        return
    qs = [q for q, _ in ctrls]
    state = "".join(str(v) for _, v in reversed(ctrls))
    qc.append(gate.control(len(qs), ctrl_state=state), qs + [target])


def _prep_rec(qc, qs, L, ctrls):
    m = len(qs)
    if m == 0 or L == 1:
        return
    if L == 2 ** m:
        for q in qs:
            _ctrl_apply(qc, HGate(), q, ctrls)
        return
    half = 2 ** (m - 1)
    if L <= half:
        _prep_rec(qc, qs[:-1], L, ctrls)
        return
    top = qs[-1]
    theta = 2 * math.acos(math.sqrt(half / L))
    _ctrl_apply(qc, RYGate(theta), top, ctrls)
    for q in qs[:-1]:
        _ctrl_apply(qc, HGate(), q, ctrls + [(top, 0)])
    _prep_rec(qc, qs[:-1], L - half, ctrls + [(top, 1)])


def uniL_recursive(n, L):
    qc = QuantumCircuit(n)
    _prep_rec(qc, list(range(n)), L, [])
    return qc


def uniL_lib(n, L):
    v = np.zeros(2 ** n); v[:L] = 1 / math.sqrt(L)
    qc = QuantumCircuit(n)
    qc.append(StatePreparation(v), range(n))
    return qc


def _uniL_spec(n, L):
    v = np.zeros(2 ** n); v[:L] = 1
    return state_spec(v)


Q05 = Task(
    id="Q05_uniform_below_L", title="uniform superposition of the first L basis states",
    source="QPC001_C1 / C2", equiv="E5", ptype="P5", params=["n", "L"],
    spec=_uniL_spec, grid=_flipL_grid,
    styles=[
        Style("recursive", uniL_recursive, "M4*", ["if-data", "mcx", "controls-grow"],
              ""),
        Style("lib", uniL_lib, "lib", ["lib"], "StatePreparation"),
    ],
    note="",
)


# =====================================================================
# Q06. Reflection 2|0><0| - I
# =====================================================================

def refl_mcx(n):
    qc = QuantumCircuit(n)
    qc.x(range(n))
    qc.h(n - 1); qc.mcx(list(range(n - 1)), n - 1); qc.h(n - 1)
    qc.x(range(n))
    return qc


def refl_ctrlstate(n):
    """MCX with ctrl_state = 0...0 instead of surrounding X gates; the Z phase by H on the target."""
    qc = QuantumCircuit(n)
    qc.x(n - 1); qc.h(n - 1)
    mcx(qc, range(n - 1), n - 1, ctrl_state="0" * (n - 1))
    qc.h(n - 1); qc.x(n - 1)
    return qc


def refl_vchain(n):
    """A Toffoli ladder into ancillas: n-2 ancillas (qubits n..2n-3), gates of bounded arity."""
    if n == 2:
        qc = QuantumCircuit(2)
        qc.x([0, 1]); qc.cz(0, 1); qc.x([0, 1])
        return qc
    anc = list(range(n, 2 * n - 2))
    qc = QuantumCircuit(2 * n - 2)
    qc.x(range(n))
    chain = []
    qc.ccx(0, 1, anc[0]); chain.append((0, 1, anc[0]))
    for i in range(2, n - 1):
        qc.ccx(i, anc[i - 2], anc[i - 1]); chain.append((i, anc[i - 2], anc[i - 1]))
    qc.cz(anc[n - 3], n - 1)
    for c1, c2, t in reversed(chain):
        qc.ccx(c1, c2, t)
    qc.x(range(n))
    return qc


def _refl_spec(n):
    M = -np.eye(2 ** n, dtype=complex); M[0, 0] = 1
    return unitary_spec(M)


Q06 = Task(
    id="Q06_reflect_zero", title="reflection 2|0><0| - I", source="QPC003_B4",
    equiv="E2 (vchain: E3)", ptype="P2", params=["n"], spec=_refl_spec,
    grid=lambda: [{"n": n} for n in range(2, 7)],
    styles=[
        Style("mcx", refl_mcx, "M2+mcx", ["mcx"], ""),
        Style("ctrlstate", refl_ctrlstate, "M0+mcx", ["mcx", "ctrl-state"], ""),
        Style("vchain", refl_vchain, "M2", ["ancilla", "two-ladders"],
              ""),
    ],
)


# =====================================================================
# Q07. Grover diffusion operator
# =====================================================================

def diff_conj(n):
    qc = QuantumCircuit(n)
    qc.h(range(n))
    qc.compose(refl_mcx(n), inplace=True)
    qc.h(range(n))
    return qc


def diff_mcp(n):
    qc = QuantumCircuit(n)
    qc.h(range(n)); qc.x(range(n))
    qc.mcp(math.pi, list(range(n - 1)), n - 1)
    qc.x(range(n)); qc.h(range(n))
    return qc


def _diff_spec(n):
    s = np.ones(2 ** n) / math.sqrt(2 ** n)
    return unitary_spec(2 * np.outer(s, s) - np.eye(2 ** n))


Q07 = Task(
    id="Q07_diffusion", title="Grover diffusion operator", source="QPC003_B5",
    equiv="E2", ptype="P2", params=["n"], spec=_diff_spec,
    grid=lambda: [{"n": n} for n in range(2, 7)],
    styles=[
        Style("conj", diff_conj, "M2+mcx", ["mcx", "compose"], ""),
        Style("mcp", diff_mcp, "M2+mcx", ["mcx"], "MCPhase(pi)"),
    ],
)


# =====================================================================
# Q08. Cyclic shift of the qubits without SwapGate
# =====================================================================

def _swap_cx(qc, a, b):
    qc.cx(a, b); qc.cx(b, a); qc.cx(a, b)


def shift_chain(n):
    qc = QuantumCircuit(n)
    for i in range(n - 1, 0, -1):
        _swap_cx(qc, i, i - 1)
    return qc


def shift_reversals(n):
    """Shift = reversal of all qubits, then reversal of qubits 1..n-1. Depth 6."""
    qc = QuantumCircuit(n)
    for i in range(n // 2):
        _swap_cx(qc, i, n - 1 - i)
    m = n - 1
    for i in range(m // 2):
        _swap_cx(qc, 1 + i, n - 1 - i)
    return qc


def _shift_spec(n):
    def f(x):
        b = bits_of(x, n)
        nb = [b[n - 1]] + b[:n - 1]
        return sum(v << i for i, v in enumerate(nb))
    return unitary_spec(perm_matrix(n, f))


Q08 = Task(
    id="Q08_cyclic_shift", title="cyclic shift of the qubits without SwapGate", source="QPC004_A3",
    equiv="E2", ptype="P2", params=["n"], spec=_shift_spec,
    grid=lambda: [{"n": n} for n in range(2, 8)],
    styles=[
        Style("chain", shift_chain, "M2", [], ""),
        Style("reversals", shift_reversals, "M3", ["fold-index"], ""),
    ],
    depth_limit=lambda n: 30, qcoder_range=lambda: [{"n": n} for n in range(2, 11)],
)


# =====================================================================
# Q09. Decrement modulo 2^n
# =====================================================================

def dec_ctrlstate(n):
    qc = QuantumCircuit(n)
    for i in range(n - 1, 0, -1):
        mcx(qc, range(i), i, ctrl_state="0" * i)
    qc.x(0)
    return qc


def _inc(qc, qs, extra_ctrls=()):
    qs = list(qs)
    for i in range(len(qs) - 1, -1, -1):
        mcx(qc, list(extra_ctrls) + qs[:i], qs[i])


def dec_xsandwich(n):
    qc = QuantumCircuit(n)
    qc.x(range(n))
    _inc(qc, range(n))
    qc.x(range(n))
    return qc


def dec_qft(n):
    qc = QuantumCircuit(n)
    qft_textbook(qc, range(n))
    for j in range(n):
        qc.p(-2 * math.pi * 2 ** j / 2 ** n, j)
    iqft_textbook(qc, range(n))
    return qc


Q09 = Task(
    id="Q09_decrement", title="decrement modulo 2^n", source="QPC004_A5",
    equiv="E2", ptype="P2", params=["n"],
    spec=lambda n: unitary_spec(perm_matrix(n, lambda x: (x - 1) % 2 ** n)),
    grid=lambda: [{"n": n} for n in range(2, 7)],
    styles=[
        Style("ctrlstate", dec_ctrlstate, "M3+mcx", ["mcx", "ctrl-state"], ""),
        Style("xsandwich", dec_xsandwich, "M3+mcx", ["mcx"], ""),
        Style("qft", dec_qft, "M3", ["angle-n"], ""),
    ],
    depth_limit=lambda n: 10, qcoder_range=lambda: [{"n": n} for n in range(2, 11)],
)


# =====================================================================
# Q10. Phase oracle exp(2 pi i c a k / 2^n)
# =====================================================================

def phase_loop(n, a):
    qc = QuantumCircuit(n + 1)
    for j in range(n):
        qc.cp(2 * math.pi * a * 2 ** j / 2 ** n, n, j)
    return qc


def phase_repeat(n, a):
    qc = QuantumCircuit(n + 1)
    for _ in range(a):
        for j in range(n):
            qc.cp(2 * math.pi * 2 ** j / 2 ** n, n, j)
    return qc


def _phase_spec(n, a):
    N = 2 ** n
    d = [np.exp(2j * np.pi * (x >> n) * a * (x % N) / N) for x in range(2 * N)]
    return unitary_spec(np.diag(d))


Q10 = Task(
    id="Q10_phase_oracle", title="phase oracle exp(2 pi i c a k / 2^n)", source="QPC004_B2",
    equiv="E2", ptype="P5", params=["n", "a"], spec=_phase_spec,
    grid=lambda: [{"n": n, "a": a} for n in range(1, 5) for a in range(1, 11)],
    styles=[
        Style("loop", phase_loop, "M3", ["angle-data"], ""),
        Style("repeat", phase_repeat, "M1(a)", ["count-a"], ""),
    ],
)


# =====================================================================
# Q11. Controlled addition of a constant a
# =====================================================================

def cadd_draper(n, a):
    qc = QuantumCircuit(n + 1)
    k, c = list(range(n)), n
    qft_textbook(qc, k)
    for j in range(n):
        qc.cp(2 * math.pi * a * 2 ** j / 2 ** n, c, k[j])
    iqft_textbook(qc, k)
    return qc


def cadd_repeat(n, a):
    qc = QuantumCircuit(n + 1)
    for _ in range(a):
        _inc(qc, range(n), extra_ctrls=[n])
    return qc


def cadd_bits(n, a):
    """For each bit t of a equal to 1: a controlled increment of qubits t..n-1."""
    qc = QuantumCircuit(n + 1)
    for t in range(n):
        if (a >> t) & 1:
            _inc(qc, range(t, n), extra_ctrls=[n])
    return qc


def _cadd_spec(n, a):
    N = 2 ** n
    return unitary_spec(perm_matrix(n + 1, lambda x: ((x % N + (x >> n) * a) % N) + N * (x >> n)))


Q11 = Task(
    id="Q11_ctrl_add_const", title="controlled addition of a constant a", source="QPC004_B3",
    equiv="E2", ptype="P5", params=["n", "a"], spec=_cadd_spec,
    grid=lambda: [{"n": n, "a": a} for n in range(1, 5) for a in range(2 ** n)],
    styles=[
        Style("draper", cadd_draper, "M3", ["angle-data"], ""),
        Style("repeat", cadd_repeat, "nonM", ["mcx", "count-a"], ""),
        Style("bits", cadd_bits, "M3*+mcx", ["mcx", "if-data"], ""),
    ],
)


# =====================================================================
# Q12. Negation modulo 2^n
# =====================================================================

def neg_xinc(n):
    qc = QuantumCircuit(n)
    qc.x(range(n))
    _inc(qc, range(n))
    return qc


def neg_or_cascade(n):
    """Bit j of -k equals k_j xor (k_0 or ... or k_(j-1))."""
    qc = QuantumCircuit(n)
    for j in range(n - 1, 0, -1):
        qc.x(j)
        mcx(qc, range(j), j, ctrl_state="0" * j)
    return qc


Q12 = Task(
    id="Q12_negate", title="negation modulo 2^n", source="QPC005_A5",
    equiv="E2", ptype="P2", params=["n"],
    spec=lambda n: unitary_spec(perm_matrix(n, lambda x: (-x) % 2 ** n)),
    grid=lambda: [{"n": n} for n in range(1, 7)],
    styles=[
        Style("xinc", neg_xinc, "M3+mcx", ["mcx"], ""),
        Style("or_cascade", neg_or_cascade, "M3+mcx", ["mcx", "ctrl-state"], ""),
    ],
    depth_limit=lambda n: 30, qcoder_range=lambda: [{"n": n} for n in range(1, 11)],
)


# =====================================================================
# Q13. m xor [k != 0]
# =====================================================================

def nz_ctrlstate(n):
    qc = QuantumCircuit(n + 1)
    m, k = 0, list(range(1, n + 1))
    mcx(qc, k, m, ctrl_state="0" * n)
    qc.x(m)
    return qc


def nz_xsandwich(n):
    qc = QuantumCircuit(n + 1)
    m, k = 0, list(range(1, n + 1))
    qc.x(k)
    qc.mcx(k, m)
    qc.x(k)
    qc.x(m)
    return qc


Q13 = Task(
    id="Q13_nonzero_flag", title="m xor [k != 0]", source="QPC005_A3",
    equiv="E2", ptype="P2", params=["n"],
    spec=lambda n: unitary_spec(perm_matrix(n + 1, lambda x: x ^ (1 if (x >> 1) != 0 else 0))),
    grid=lambda: [{"n": n} for n in range(1, 7)],
    styles=[
        Style("ctrlstate", nz_ctrlstate, "M0+mcx", ["mcx", "ctrl-state"], ""),
        Style("xsandwich", nz_xsandwich, "M2+mcx", ["mcx"], ""),
    ],
    depth_limit=lambda n: 4, qcoder_range=lambda: [{"n": n} for n in range(1, 11)],
)


# =====================================================================
# Q14. X error on qubit k-1
# =====================================================================

def err_ifchain(k):
    qc = QuantumCircuit(3)
    if k == 1:
        qc.x(0)
    elif k == 2:
        qc.x(1)
    elif k == 3:
        qc.x(2)
    return qc


def err_index(k):
    qc = QuantumCircuit(3)
    if k > 0:
        qc.x(k - 1)
    return qc


def _err_spec(k):
    return unitary_spec(perm_matrix(3, lambda x: x ^ (1 << (k - 1)) if k else x))


Q14 = Task(
    id="Q14_error_k", title="X error on qubit k-1", source="QPC006_A1",
    equiv="E2", ptype="P5 (data without n)", params=["k"], spec=_err_spec,
    grid=lambda: [{"k": k} for k in range(4)],
    styles=[
        Style("ifchain", err_ifchain, "M0", ["if-param"], ""),
        Style("index", err_index, "M0", ["if-param"], ""),
    ],
    note="",
)


TASKS = [Q01, Q02, Q03, Q04, Q05, Q06, Q07, Q08, Q09, Q10, Q11, Q12, Q13, Q14]
