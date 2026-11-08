"""Textbook tasks of the own generator bank: GHZ, W, BV, DJ, QFT, inverse QFT, Grover, QPE, adder, repetition
code, and two tasks whose parameter is a repetition count (P1)."""
from __future__ import annotations

import math
import random

import numpy as np
from qiskit import QuantumCircuit
from qiskit.circuit.library import QFTGate, StatePreparation

from common import (Style, Task, bits_of, dft_matrix, iqft_textbook, mcx,
                    qft_textbook, state_spec, unitary_spec)


# =====================================================================
# C01. GHZ
# =====================================================================

def ghz_ladder(n):
    qc = QuantumCircuit(n)
    qc.h(0)
    for i in range(n - 1):
        qc.cx(i, i + 1)
    return qc


def ghz_star(n):
    qc = QuantumCircuit(n)
    qc.h(0)
    for i in range(1, n):
        qc.cx(0, i)
    return qc


def ghz_tree(n):
    """Doubling: after step t the qubits 0..have-1 are entangled and are copied to have..2*have-1."""
    qc = QuantumCircuit(n)
    qc.h(0)
    have = 1
    while have < n:
        step = min(have, n - have)
        for i in range(step):
            qc.cx(i, have + i)
        have += step
    return qc


def ghz_spec(n):
    v = np.zeros(2 ** n); v[0] = v[-1] = 1
    return state_spec(v)


C01 = Task(
    id="C01_ghz", title="GHZ state", source="classic (cf. QPC002_A3)",
    equiv="E5", ptype="P2", params=["n"], spec=ghz_spec,
    grid=lambda: [{"n": n} for n in range(2, 8)],
    styles=[
        Style("ladder", ghz_ladder, "M2", [], ""),
        Style("star", ghz_star, "M3", [], ""),
        Style("tree", ghz_tree, "M4-tree", ["log-depth", "index-2^t"],
              ""),
    ],
)


# =====================================================================
# C02. W state
# =====================================================================

def w_cascade(n):
    qc = QuantumCircuit(n)
    qc.x(0)
    for i in range(n - 1):
        theta = 2 * math.acos(math.sqrt(1 / (n - i)))
        qc.cry(theta, i, i + 1)
        qc.cx(i + 1, i)
    return qc


def _w_split(qc, lo, size):
    """The excitation is on qubit lo; spread it evenly over lo..lo+size-1."""
    if size <= 1:
        return
    a = size // 2
    b = size - a
    theta = 2 * math.acos(math.sqrt(a / size))
    qc.cry(theta, lo, lo + a)
    qc.cx(lo + a, lo)
    _w_split(qc, lo, a)
    _w_split(qc, lo + a, b)


def w_tree(n):
    qc = QuantumCircuit(n)
    qc.x(0)
    _w_split(qc, 0, n)
    return qc


def w_lib(n):
    v = np.zeros(2 ** n)
    for i in range(n):
        v[1 << i] = 1 / math.sqrt(n)
    qc = QuantumCircuit(n)
    qc.append(StatePreparation(v), range(n))
    return qc


def w_spec(n):
    v = np.zeros(2 ** n)
    for i in range(n):
        v[1 << i] = 1
    return state_spec(v)


C02 = Task(
    id="C02_w", title="W state", source="classic = QPC003_A4",
    equiv="E5", ptype="P2", params=["n"], spec=w_spec,
    grid=lambda: [{"n": n} for n in range(2, 8)],
    styles=[
        Style("cascade", w_cascade, "M2~", ["angle-n"],
              ""),
        Style("tree", w_tree, "M4-tree", ["angle-n", "log-depth", "two-calls"],
              ""),
        Style("lib", w_lib, "lib", ["lib"], ""),
    ],
    depth_limit=lambda n: 10, qcoder_range=lambda: [{"n": n} for n in range(2, 16)],
    note="",
)


# =====================================================================
# C03. Bernstein-Vazirani (whole algorithm)
# =====================================================================

def bv_phase(n, s):
    qc = QuantumCircuit(n)
    qc.h(range(n))
    for i in range(n):
        if (s >> i) & 1:
            qc.z(i)
    qc.h(range(n))
    return qc


def bv_ancilla(n, s):
    qc = QuantumCircuit(n + 1)
    y = n
    qc.x(y); qc.h(y)
    qc.h(range(n))
    for i in range(n):
        if (s >> i) & 1:
            qc.cx(i, y)
    qc.h(range(n))
    qc.h(y); qc.x(y)
    return qc


def bv_shortcut(n, s):
    """The circuit of the whole algorithm equals X^s, since H Z H = X: a shortcut solution."""
    qc = QuantumCircuit(n)
    for i in range(n):
        if (s >> i) & 1:
            qc.x(i)
    return qc


def bv_spec(n, s):
    N = 2 ** n
    M = np.zeros((N, N), dtype=complex)
    for x in range(N):
        M[x ^ s, x] = 1
    return unitary_spec(M)


def _bv_grid():
    out = []
    for n in range(1, 6):
        ss = range(2 ** n) if n <= 3 else random.Random(n).sample(range(2 ** n), 6)
        out += [{"n": n, "s": s} for s in ss]
    return out


C03 = Task(
    id="C03_bv", title="Bernstein-Vazirani (whole algorithm)", source="classic",
    equiv="E2 (style ancilla: E3)", ptype="P5", params=["n", "s"], spec=bv_spec, grid=_bv_grid,
    styles=[
        Style("phase", bv_phase, "M3*", ["if-data"], ""),
        Style("ancilla", bv_ancilla, "M3*", ["if-data", "ancilla"],
              ""),
        Style("shortcut", bv_shortcut, "M3*", ["if-data"], ""),
    ],
    note="",
)


# =====================================================================
# C04. Deutsch-Jozsa with the oracle as an argument
# =====================================================================

def dj_oracle(n, table):
    """O_f |x>|y> = |x>|y xor f(x)>, y is qubit n."""
    o = QuantumCircuit(n + 1)
    for x, fx in enumerate(table):
        if fx:
            mcx(o, range(n), n, ctrl_state="".join(str(b) for b in bits_of(x, n)))
    return o


def dj_kickback(n, o):
    qc = QuantumCircuit(n + 1)
    qc.x(n); qc.h(n)
    qc.h(range(n))
    qc.compose(o, range(n + 1), inplace=True)
    qc.h(range(n))
    qc.h(n); qc.x(n)
    return qc


def dj_compute_uncompute(n, o):
    """Two calls of the oracle: f into an ancilla, Z, the oracle again."""
    qc = QuantumCircuit(n + 1)
    qc.h(range(n))
    qc.compose(o, range(n + 1), inplace=True)
    qc.z(n)
    qc.compose(o, range(n + 1), inplace=True)
    qc.h(range(n))
    return qc


def dj_spec(n, o, _table):
    N = 2 ** n
    Hn = np.array([[1]], dtype=complex)
    H = np.array([[1, 1], [1, -1]]) / math.sqrt(2)
    for _ in range(n):
        Hn = np.kron(H, Hn)
    D = np.diag([(-1) ** f for f in _table]).astype(complex)
    M = Hn @ D @ Hn
    return M, n, n


def _dj_grid():
    out = []
    for n in range(1, 5):
        N = 2 ** n
        rng = random.Random(100 + n)
        tables = [[0] * N, [1] * N]
        for _ in range(3):
            t = [1] * (N // 2) + [0] * (N - N // 2)
            rng.shuffle(t)
            tables.append(t)
        for t in tables:
            out.append({"n": n, "o": dj_oracle(n, t), "_table": t})
    return out


C04 = Task(
    id="C04_dj", title="Deutsch-Jozsa with the oracle as an argument", source="classic (cf. QPC003_B2)",
    equiv="E3", ptype="outside P1-P5: a circuit argument", params=["n", "o"], spec=dj_spec, grid=_dj_grid,
    styles=[
        Style("kickback", dj_kickback, "HO", ["oracle-arg", "ancilla"], ""),
        Style("compute_uncompute", dj_compute_uncompute, "HO", ["oracle-arg", "ancilla"],
              ""),
    ],
    note="",
)


# =====================================================================
# C05. Quantum Fourier transform
# =====================================================================

def qft_loops(n):
    qc = QuantumCircuit(n)
    qft_textbook(qc, range(n))
    return qc


def _qft_rec(qc, qs):
    if not qs:
        return
    top = qs[-1]
    qc.h(top)
    for k in reversed(range(len(qs) - 1)):
        qc.cp(math.pi / 2 ** (len(qs) - 1 - k), top, qs[k])
    _qft_rec(qc, qs[:-1])


def qft_recursive(n):
    qc = QuantumCircuit(n)
    _qft_rec(qc, list(range(n)))
    for i in range(n // 2):
        a, b = i, n - 1 - i
        qc.cx(a, b); qc.cx(b, a); qc.cx(a, b)
    return qc


def qft_lib(n):
    qc = QuantumCircuit(n)
    qc.append(QFTGate(n), range(n))
    return qc


C05 = Task(
    id="C05_qft", title="quantum Fourier transform", source="classic = QPC002_B4",
    equiv="E2", ptype="P2", params=["n"], spec=lambda n: unitary_spec(dft_matrix(n, +1)),
    grid=lambda: [{"n": n} for n in range(1, 7)],
    styles=[
        Style("loops", qft_loops, "M3", ["angle-family"], ""),
        Style("recursive", qft_recursive, "M4", ["angle-family"],
              ""),
        Style("lib", qft_lib, "lib", ["lib"], ""),
    ],
    depth_limit=lambda n: 25, qcoder_range=lambda: [{"n": n} for n in range(1, 11)],
    note="",
)


def iqft_loops(n):
    qc = QuantumCircuit(n)
    iqft_textbook(qc, range(n))
    return qc


def iqft_inverse(n):
    return qft_loops(n).inverse()


C06 = Task(
    id="C06_iqft", title="inverse QFT", source="classic",
    equiv="E2", ptype="P2", params=["n"], spec=lambda n: unitary_spec(dft_matrix(n, -1)),
    grid=lambda: [{"n": n} for n in range(1, 7)],
    styles=[
        Style("loops", iqft_loops, "M3", ["angle-family"], ""),
        Style("inverse", iqft_inverse, "lib", ["inverse()"], ""),
    ],
)


# =====================================================================
# C07. Grover search (n, w)
# =====================================================================

def grover_iters(n):
    return int(math.floor(math.pi / 4 * math.sqrt(2 ** n)))


def _oracle_mcx(qc, n, w):
    for i in range(n):
        if not (w >> i) & 1:
            qc.x(i)
    if n == 1:
        qc.z(0)
    else:
        qc.h(n - 1); qc.mcx(list(range(n - 1)), n - 1); qc.h(n - 1)
    for i in range(n):
        if not (w >> i) & 1:
            qc.x(i)


def _diffusion_mcx(qc, n):
    qc.h(range(n)); qc.x(range(n))
    if n == 1:
        qc.z(0)
    else:
        qc.h(n - 1); qc.mcx(list(range(n - 1)), n - 1); qc.h(n - 1)
    qc.x(range(n)); qc.h(range(n))


def grover_mcx(n, w):
    qc = QuantumCircuit(n)
    qc.h(range(n))
    for _ in range(grover_iters(n)):
        _oracle_mcx(qc, n, w)
        _diffusion_mcx(qc, n)
    return qc


def grover_mcp(n, w):
    """Oracle and diffusion by a multi-controlled phase with ctrl_state."""
    qc = QuantumCircuit(n)
    qc.h(range(n))
    for _ in range(grover_iters(n)):
        if n == 1:
            qc.p(math.pi, 0) if w == 1 else (qc.x(0), qc.p(math.pi, 0), qc.x(0))
        else:
            for i in range(n):
                if not (w >> i) & 1:
                    qc.x(i)
            qc.mcp(math.pi, list(range(n - 1)), n - 1)
            for i in range(n):
                if not (w >> i) & 1:
                    qc.x(i)
        qc.h(range(n)); qc.x(range(n))
        if n == 1:
            qc.p(math.pi, 0)
        else:
            qc.mcp(math.pi, list(range(n - 1)), n - 1)
        qc.x(range(n)); qc.h(range(n))
    return qc


def _grover_matrix(n, w, iters):
    N = 2 ** n
    s = np.ones(N) / math.sqrt(N)
    O = np.eye(N); O[w, w] = -1
    D = 2 * np.outer(s, s) - np.eye(N)
    Hn = np.array([[1]], dtype=complex)
    H = np.array([[1, 1], [1, -1]]) / math.sqrt(2)
    for _ in range(n):
        Hn = np.kron(H, Hn)
    M = Hn.astype(complex)
    for _ in range(iters):
        M = D @ O @ M
    return M


def _grover_grid():
    out = []
    for n in range(2, 6):
        ws = range(2 ** n) if n <= 3 else random.Random(n).sample(range(2 ** n), 4)
        out += [{"n": n, "w": w} for w in ws]
    return out


C07 = Task(
    id="C07_grover", title="Grover search (n, w)", source="classic",
    equiv="E2", ptype="P5", params=["n", "w"],
    spec=lambda n, w: unitary_spec(_grover_matrix(n, w, grover_iters(n))), grid=_grover_grid,
    styles=[
        Style("mcx", grover_mcx, "nonM", ["mcx", "count-sqrt", "if-data"],
              ""),
        Style("mcp", grover_mcp, "nonM", ["mcx", "count-sqrt", "if-data"], ""),
    ],
    note="",
)


# =====================================================================
# C08. Grover at n = 3, w = 5; the parameter is the number of iterations k
# =====================================================================

G_N, G_W = 3, 5


def _grover_step():
    st = QuantumCircuit(G_N)
    _oracle_mcx(st, G_N, G_W)
    _diffusion_mcx(st, G_N)
    return st


def groverk_loop(k):
    qc = QuantumCircuit(G_N)
    qc.h(range(G_N))
    st = _grover_step()
    for _ in range(k):
        qc.compose(st, inplace=True)
    return qc


def groverk_power(k):
    qc = QuantumCircuit(G_N)
    qc.h(range(G_N))
    qc.compose(_grover_step().repeat(k).decompose(), inplace=True) if k > 0 else None
    return qc


C08 = Task(
    id="C08_grover_k", title="Grover at n = 3, w = 5; the parameter is the number of iterations k", source="classic",
    equiv="E2", ptype="P1", params=["k"],
    spec=lambda k: unitary_spec(_grover_matrix(G_N, G_W, k)),
    grid=lambda: [{"k": k} for k in range(0, 9)],
    styles=[
        Style("loop", groverk_loop, "M1", [], ""),
        Style("repeat", groverk_power, "lib", ["repeat()"], "QuantumCircuit.repeat(k)"),
    ],
)


# =====================================================================
# C09. Trotterized Ising chain, n = 4; the parameter is the number of steps k
# =====================================================================

T_N, T_ZZ, T_X = 4, 0.3, 0.2


def _zz(qc, i):
    qc.rzz(2 * T_ZZ, i, i + 1)


def trotter_seq(k):
    qc = QuantumCircuit(T_N)
    for _ in range(k):
        for i in range(T_N - 1):
            _zz(qc, i)
        qc.rx(2 * T_X, range(T_N))
    return qc


def trotter_brick(k):
    """The ZZ terms commute: first the even pairs, then the odd pairs."""
    qc = QuantumCircuit(T_N)
    for _ in range(k):
        for i in range(0, T_N - 1, 2):
            _zz(qc, i)
        for i in range(1, T_N - 1, 2):
            _zz(qc, i)
        qc.rx(2 * T_X, range(T_N))
    return qc


def trotter_pairs(k):
    """A loop over pairs of steps plus a remainder for odd k: the structure depends on k mod 2."""
    def step(qc):
        for i in range(T_N - 1):
            _zz(qc, i)
        qc.rx(2 * T_X, range(T_N))
    qc = QuantumCircuit(T_N)
    for _ in range(k // 2):
        step(qc); step(qc)
    if k % 2 == 1:
        step(qc)
    return qc


def _trotter_matrix(k):
    N = 2 ** T_N
    diag = np.zeros(N, dtype=complex)
    for x in range(N):
        z = [1 - 2 * ((x >> i) & 1) for i in range(T_N)]
        diag[x] = np.exp(-1j * T_ZZ * sum(z[i] * z[i + 1] for i in range(T_N - 1)))
    ZZ = np.diag(diag)
    rx = np.array([[math.cos(T_X), -1j * math.sin(T_X)], [-1j * math.sin(T_X), math.cos(T_X)]])
    RX = np.array([[1]], dtype=complex)
    for _ in range(T_N):
        RX = np.kron(rx, RX)
    step = RX @ ZZ
    return np.linalg.matrix_power(step, k)


C09 = Task(
    id="C09_trotter_k", title="Trotterized Ising chain, n = 4; the parameter is the number of steps k",
    source="classic", equiv="E2", ptype="P1 (style pairs: P4)", params=["k"],
    spec=lambda k: unitary_spec(_trotter_matrix(k)), grid=lambda: [{"k": k} for k in range(0, 7)],
    styles=[
        Style("seq", trotter_seq, "M1", ["angle-const"], ""),
        Style("brick", trotter_brick, "M1", ["angle-const"], ""),
        Style("pairs", trotter_pairs, "M1+mod2", ["if-parity"], ""),
    ],
)


# =====================================================================
# C10. Phase estimation of P(2 pi * 0.3); the parameter is the number of counting qubits t
# =====================================================================

PHI = 0.3


def qpe_angles(t):
    qc = QuantumCircuit(t + 1)
    qc.x(t)
    qc.h(range(t))
    for j in range(t):
        qc.cp(2 * math.pi * PHI * 2 ** j, j, t)
    iqft_textbook(qc, range(t))
    return qc


def qpe_repeat(t):
    """Controlled U^(2^j) as 2^j copies of CP(2 pi phi): the gate count grows as 2^t."""
    qc = QuantumCircuit(t + 1)
    qc.x(t)
    qc.h(range(t))
    for j in range(t):
        for _ in range(2 ** j):
            qc.cp(2 * math.pi * PHI, j, t)
    qc.compose(qft_loops(t).inverse(), range(t), inplace=True)
    return qc


def qpe_spec(t):
    N = 2 ** t
    y = np.arange(N)
    v = np.exp(2j * np.pi * PHI * y) / math.sqrt(N)
    out = dft_matrix(t, -1) @ v
    full = np.zeros(2 * N, dtype=complex)
    full[N:] = out
    return state_spec(full)


C10 = Task(
    id="C10_qpe", title="phase estimation of P(2 pi * 0.3); the parameter is the number of counting qubits t",
    source="classic", equiv="E5", ptype="P2", params=["t"], spec=qpe_spec,
    grid=lambda: [{"t": t} for t in range(1, 6)],
    styles=[
        Style("angles", qpe_angles, "M3", ["angle-n"],
              ""),
        Style("repeat", qpe_repeat, "nonM", ["count-2^j", "inverse()"],
              ""),
    ],
)


# =====================================================================
# C11. Adder a + b mod 2^n
# =====================================================================

def _maj(qc, x, y, z):
    qc.cx(z, y); qc.cx(z, x); qc.ccx(x, y, z)


def _uma(qc, x, y, z):
    qc.ccx(x, y, z); qc.cx(z, x); qc.cx(x, y)


def adder_cuccaro(n):
    """Cuccaro et al.: a ladder of MAJ, then the reverse ladder of UMA; one ancilla (qubit 2n)."""
    qc = QuantumCircuit(2 * n + 1)
    a = list(range(n)); b = list(range(n, 2 * n)); c = 2 * n
    _maj(qc, c, b[0], a[0])
    for i in range(1, n):
        _maj(qc, a[i - 1], b[i], a[i])
    for i in reversed(range(1, n)):
        _uma(qc, a[i - 1], b[i], a[i])
    _uma(qc, c, b[0], a[0])
    return qc


def adder_draper(n):
    """QFT adder: QFT on b, phases 2 pi 2^(j+k)/2^n controlled by a_j, inverse QFT."""
    qc = QuantumCircuit(2 * n)
    a = list(range(n)); b = list(range(n, 2 * n))
    qft_textbook(qc, b)
    for j in range(n):
        for k in range(n - j):
            qc.cp(2 * math.pi * 2 ** (j + k) / 2 ** n, a[j], b[k])
    iqft_textbook(qc, b)
    return qc


def adder_spec(n):
    N = 2 ** n
    M = np.zeros((N * N, N * N), dtype=complex)
    for x in range(N * N):
        av, bv = x % N, x // N
        M[av + N * ((av + bv) % N), x] = 1
    return unitary_spec(M)


C11 = Task(
    id="C11_adder", title="adder a + b mod 2^n", source="classic",
    equiv="E2 (cuccaro: E3)", ptype="P2", params=["n"], spec=adder_spec,
    grid=lambda: [{"n": n} for n in range(1, 5)],
    styles=[
        Style("cuccaro", adder_cuccaro, "M2", ["ancilla", "two-ladders"],
              ""),
        Style("draper", adder_draper, "M3", ["angle-family"], ""),
    ],
    note="",
)


# =====================================================================
# C12. Encoder of the repetition code
# =====================================================================

def rep_ladder(n):
    qc = QuantumCircuit(n)
    for i in range(n - 1):
        qc.cx(i, i + 1)
    return qc


def rep_star(n):
    qc = QuantumCircuit(n)
    for i in range(1, n):
        qc.cx(0, i)
    return qc


def rep_tree(n):
    qc = QuantumCircuit(n)
    have = 1
    while have < n:
        step = min(have, n - have)
        for i in range(step):
            qc.cx(i, have + i)
        have += step
    return qc


def rep_spec(n):
    W = np.zeros((2 ** n, 2), dtype=complex)
    W[0, 0] = 1; W[2 ** n - 1, 1] = 1
    return W, 1, n


C12 = Task(
    id="C12_repetition", title="encoder of the repetition code", source="classic",
    equiv="E3", ptype="P2", params=["n"], spec=rep_spec,
    grid=lambda: [{"n": n} for n in range(2, 8)],
    styles=[
        Style("ladder", rep_ladder, "M2", [], "CX(i, i+1)"),
        Style("star", rep_star, "M3", [], "CX(0, i)"),
        Style("tree", rep_tree, "M4-tree", ["log-depth", "index-2^t"], ""),
    ],
    note="",
)


TASKS = [C01, C02, C03, C04, C05, C06, C07, C08, C09, C10, C11, C12]
