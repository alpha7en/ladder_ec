"""Ladder programs (block ladders) and their automata (Sec. 4 and Thm. 1 of the paper).

A program is a list of segments on n + e qubits (e ancillas at the end):
  Fixed(gates)                         gates on qubits k or n + k (a fixed segment),
  Loop(start, stop, step, body)        for i in range(start, stop, step): body, any step != 0,
                                       body gates on qubits i + o (the window) and on fixed qubits.
A reference is ("c", k) = k, ("n", k) = n + k or ("i", k) = i + k. A gate may carry conditions
((ref, v), ...): it is applied only if a[ref] = v for each of them (nested `if` on bits).

Automaton of a loop (step +1). Block B_i acts on the window i, ..., i + w - 1 and on the fixed qubits. Reading the
qubits left to right, B_i is applied when the letter of qubit i + w - 1 is read: then all its inputs are known.
The state keeps the current values of the open window qubits and of the fixed qubits (one basis vector each, the
"carry"), the x-letters of the open qubits (to project qubit i onto <x_i| once B_i is done) and the bits a
the conditions need. Fixed qubits at the end of the register (the set E of Definition 3) are read after the loop,
so their input y is guessed at the start and checked when it is read. With step s, the first min(s, w) window
qubits are finished after each block, and a counter tells when the next block is due. Hence D <= (s + 1) 4^{w-1}
4^{h} up to the bits and the two ends (Thm. 1). A loop with a negative step is the mirror image of one with a
positive step, and its automaton is the transposed one (the word is read backwards). Proof: Appendix A.
"""
from collections import namedtuple
from functools import lru_cache
import itertools
import numpy as np
from .field import apply
from .automata import explore, minimize, compose, WA

Gate = namedtuple("Gate", "op qubits angle cond", defaults=(None, ()))
Fixed = namedtuple("Fixed", "gates")
Loop = namedtuple("Loop", "start stop step body")
Family = namedtuple("Family", "segments extra n_min")


def c(k): return ("c", k)
def n(k=0): return ("n", k)
def i(k=0): return ("i", k)


def letter(l):
    return l >> 2, (l >> 1) & 1, l & 1                     # (a, x, y)


# ------------------------------------------------------------------ direct semantics (for small n)
def pos(ref, n_, i_=None):
    kind, k = ref
    return {"c": k, "n": n_ + k, "i": (i_ or 0) + k}[kind]


def expand(fam, n_, a):
    """the gate list (op, qubits, angle) of the circuit for size n and bits a (a list of length >= n + e)"""
    out = []
    def emit(g, i_=None):
        if all(a[pos(r, n_, i_)] == v for r, v in g.cond):
            out.append((g.op, [pos(q, n_, i_) for q in g.qubits], g.angle))
    for seg in fam.segments:
        if isinstance(seg, Fixed):
            for g in seg.gates:
                emit(g)
        else:
            for i_ in range(pos(seg.start, n_), pos(seg.stop, n_), seg.step):
                for g in seg.body:
                    emit(g, i_)
    return out


def matrix(F, gates, nq):
    """the circuit matrix as an array of shape (2,)*nq + (2,)*nq: axes = output qubits, then input qubits"""
    psi = np.eye(2 ** nq, dtype=np.int64).reshape((2,) * nq + (2 ** nq,))
    for op, qs, t in gates:
        psi = apply(F, F.gate(op, t), qs, psi)
    return psi.reshape((2,) * 2 * nq)


# ------------------------------------------------------------------ automata of segments
def fixed_wa(F, seg, e):
    """Fixed gates touch positions k (start) and N + d (end, d < 0). The machine collects the letters of the start
    positions and guesses those of the end positions, weighs the word by the entry of the gate product as soon as
    the last start letter is read, carries only the guess through the gap of identity letters and checks it at the
    end (as for the end fixed qubits of a loop). So the gap holds 8^|end| states, not 8^(|start| + |end|)."""
    refs = sorted({q for g in seg.gates for q in g.qubits} | {r for g in seg.gates for r, _ in g.cond})
    key = lambda r: ("c", r[1]) if r[0] == "c" else ("n", r[1] - e)          # n + k = N + (k - e)
    keys = [key(r) for r in refs]
    S = [k for kind, k in keys if kind == "c"]
    E = [("n", k) for kind, k in keys if kind == "n"]
    s_len, e_len = (max(S) + 1 if S else 0), (-min(k for _, k in E) if E else 0)
    qubits = sorted({key(q) for g in seg.gates for q in g.qubits})

    @lru_cache(None)
    def U(abits):                                          # gate product for the collected bits a
        bit = dict(abits)
        gs = [(g.op, [qubits.index(key(q)) for q in g.qubits], g.angle) for g in seg.gates
              if all(bit[key(r)] == v for r, v in g.cond)]
        return matrix(F, gs, len(qubits))

    def weight(col):                                       # col: (position, a, x, y) of every touched position
        got = {k: (a, x, y) for k, a, x, y in col}
        if any(got[k][1] != got[k][2] for k in got if k not in qubits):
            return 0
        M = U(tuple(sorted((k, v[0]) for k, v in got.items())))
        return int(M[tuple(got[q][1] for q in qubits) + tuple(got[q][2] for q in qubits)])

    def gap(col):                                          # after the last start letter: weigh and guess the end
        out = []
        for guess in itertools.product(*[[(k,) + letter(l) for l in range(8)] for k in E]):
            w = weight(col + guess)
            if w:
                out.append((("G", 0, guess), w))
        return out

    def read(k, l, col):                                   # start letter l at position ("c", k)
        a, x, y = letter(l)
        if ("c", k) in keys:
            return [col + ((("c", k), a, x, y),)]
        return [col] if x == y else []

    def check(d, l, guess):                                # end letter l at position ("n", d)
        a, x, y = letter(l)
        return (("n", d), a, x, y) in guess if ("n", d) in E else x == y

    def step(s, l):
        ph, k, col = s
        if ph == "S":
            if k + 1 < s_len:
                return [(("S", k + 1, c2), 1) for c2 in read(k, l, col)]
            return [t for c2 in read(k, l, col) for t in gap(c2)]
        if ph == "G":
            out = [(s, 1)] if letter(l)[1] == letter(l)[2] else []
            if e_len and check(-e_len, l, col):
                out.append((("E", 1, col), 1))
            return out
        if k < e_len and check(k - e_len, l, col):
            return [(("E", k + 1, col), 1)]
        return []

    final = lambda s: int(s[:2] == (("E", e_len) if e_len else ("G", 0)))
    return explore(F.p, 8, [(("S", 0, ()), 1)] if s_len else gap(()), step, final)


def mirror(seg, e):
    """the loop read from the other end of the register: k <-> N - 1 - k"""
    m = lambda r: {"c": ("n", e - 1 - r[1]), "n": ("c", e - 1 - r[1]), "i": ("i", -r[1])}[r[0]]
    body = [Gate(g.op, [m(q) for q in g.qubits], g.angle, tuple((m(r), v) for r, v in g.cond)) for g in seg.body]
    return Loop(m(seg.start), m(seg.stop), -seg.step, body)


def loop_wa(F, seg, e):
    if seg.step < 0:
        W = loop_wa(F, mirror(seg, e), e)
        return WA(W.beta, W.A.transpose(0, 2, 1), W.alpha)
    assert seg.start[0] == "c" and seg.stop[0] == "n"
    st = seg.step
    refs = [q for g in seg.body for q in g.qubits] + [r for g in seg.body for r, _ in g.cond]
    offs = [k for kind, k in refs if kind == "i"]
    lo, w = min(offs), max(offs) - min(offs) + 1
    i0, dd = seg.start[1] + lo, seg.stop[1] + lo - e      # blocks i = i0, i0 + st, ... < N + dd
    K = 1 - dd - w                                          # letters after the last block: K, ..., K + st - 1
    SF = sorted({k for kind, k in refs if kind == "c"})
    EF = sorted({k - e for kind, k in refs if kind == "n"})
    assert K >= 0 and all(k < i0 for k in SF) and all(K + d >= 0 for d in EF)
    assert not any(r[0] == "n" for g in seg.body for r, _ in g.cond), "condition on a fixed bit at the end of the register"
    wire = lambda r: {"i": lambda k: k - lo, "c": lambda k: w + SF.index(k),
                      "n": lambda k: w + len(SF) + EF.index(k - e)}[r[0]](r[1])
    nq = w + len(SF) + len(EF)
    need = {wire(r) for g in seg.body for r, _ in g.cond}      # wires whose bit a some condition reads
    # a letter on window wire j moves to j - st, j - 2 st, ... with the next blocks, a start fixed qubit stays: the bit
    # is kept only if a condition of this or a later block reads it, otherwise 0 is stored in its place
    keep = lambda j, a: a if (j in need if j >= w else any(j - t in need for t in range(0, j + 1, st))) else 0

    @lru_cache(None)
    def block(abits):                                      # abits: bits of the window, then of the start fixed qubits
        gs = [(g.op, [wire(q) for q in g.qubits], g.angle) for g in seg.body
              if all(abits[wire(r)] == v for r, v in g.cond)]
        return matrix(F, gs, nq).reshape(2 ** nq, 2 ** nq)

    def step(s, l):
        a, x, y = letter(l)
        ph, k, win, sh, eh = s
        if ph == "P":                                      # before the window: start fixed qubits or identity
            if k in SF:
                sh = sh + ((x, keep(w + len(sh), a), y),)
            elif x != y:
                return []
            return [((("P", k + 1) if k + 1 < i0 else ("R", w)) + (win, sh, eh), 1)]
        if ph == "R":                                      # k = letters until the next block
            if k > w:                                      # between two windows (step > w): identity
                return [(("R", k - 1, win, sh, eh), 1)] if x == y else []
            win = win + ((x, keep(len(win), a), y),)
            if k > 1:
                return [(("R", k - 1, win, sh, eh), 1)]
            col = block(tuple(v[1] for v in win + sh))
            idx = int("".join(str(v[2]) for v in win + sh) + "".join(str(v[1]) for v in eh) or "0", 2)
            out = []
            for o in np.flatnonzero(col[:, idx]):
                bits = [int(b) for b in format(o, "0%db" % nq)]
                f = min(st, w)                             # qubits i, ..., i + f - 1 are finished:
                if any(b != v[0] for v, b in zip(win[:f], bits[:f])):   # project them onto <x|
                    continue
                win2 = tuple((v[0], keep(j, v[1]), b) for j, (v, b) in enumerate(zip(win[f:], bits[f:w])))
                sh2 = tuple((v[0], v[1], b) for v, b in zip(sh, bits[w:w + len(SF)]))
                eh2 = tuple((g, b) for (g, _), b in zip(eh, bits[w + len(SF):]))
                wt = int(col[o, idx])
                out.append((("R", st, win2, sh2, eh2), wt))
                if all(v[0] == v[2] for v in win2 + sh2):  # last block: project the rest, go to the suffix,
                    out += [(("S", (0, K + r), (), (), eh2), wt) for r in range(st)]   # whose length is K + (N + dd - 1 - i) mod st
            return out
        j, KK = k
        if j >= KK:
            return []
        d = j - KK
        if d in EF:
            g, cbit = eh[EF.index(d)]
            if (g, cbit) != (y, x):
                return []
        elif x != y:
            return []
        return [(("S", (j + 1, KK), (), (), eh), 1)]

    final = lambda s: int(s[0] == "S" and s[1][0] == s[1][1])
    init = [((("P", 0) if i0 else ("R", w)) + ((), (), tuple((g, g) for g in gs)), 1)
            for gs in itertools.product((0, 1), repeat=len(EF))]
    return explore(F.p, 8, init, step, final)


def segment_wa(F, seg, e):
    return fixed_wa(F, seg, e) if isinstance(seg, Fixed) else loop_wa(F, seg, e)


def family_wa(F, fam):
    """automaton of the whole family: product of the segment automata in circuit order, minimised"""
    W = None
    for seg in fam.segments:
        S = minimize(F.p, segment_wa(F, seg, fam.extra))
        W = S if W is None else minimize(F.p, compose(F.p, S, W))
    return W
