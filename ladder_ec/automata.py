"""Weighted automata over F_p (Sec. 5 of the paper).

An automaton W = (alpha, A, beta) with letters l = 0..k-1 computes f(w) = alpha A(w_1) ... A(w_N) beta.
A circuit family is read qubit by qubit: the letter of qubit j is (a_j, x_j, y_j), letter index 4a + 2x + y,
and f(w) = <x| U_{n,a} |y>.

Minimisation and the zero test are the forward/backward reductions of Kiefer (2020), Prop. 2.1-2.2, Thm. 3.6.
All products are int64 with one reduction mod p after the sum over the inner index; no_overflow guards that.
"""
from collections import deque, namedtuple
import numpy as np

WA = namedtuple("WA", "alpha A beta")
USED = {"max_p": 0, "max_dim": 0, "max_load": 0}          # largest prime, inner dimension and D (p - 1)^2 seen


def no_overflow(p, D):
    """int64 products below are reduced mod p only after the sum over an inner index of length D: exact iff
    D (p - 1)^2 < 2^63 (p < 2^25 allows D < 2^13). Checked before every such product (Sec. 6.2); the running
    maxima are reported in the field `used` of every output row."""
    load = D * (p - 1) ** 2
    assert load < 2 ** 63, "int64 overflow: inner dimension %d with p = %d gives D (p-1)^2 = %d >= 2^63" % (D, p, load)
    if p > USED["max_p"]:
        USED["max_p"] = p
    if D > USED["max_dim"]:
        USED["max_dim"] = D
    if load > USED["max_load"]:
        USED["max_load"] = load


def explore(p, letters, init, step, final):
    """automaton of a finite machine: init = [(state, weight)], step(state, letter) = [(state', weight)],
    final(state) = weight. States are hashable; only reachable ones are kept."""
    index, todo, edges = {}, deque(), []
    def add(s):
        if s not in index:
            index[s] = len(index)
            todo.append(s)
        return index[s]
    starts = [(add(s), w) for s, w in init]
    while todo:
        s = todo.popleft()
        for l in range(letters):
            for t, w in step(s, l):
                if w % p:
                    edges.append((l, index[s], add(t), w))
    D = len(index)
    alpha, beta, A = np.zeros(D, np.int64), np.zeros(D, np.int64), np.zeros((letters, D, D), np.int64)
    for i, w in starts:
        alpha[i] = (alpha[i] + w) % p
    for s, i in index.items():
        beta[i] = final(s) % p
    for l, i, j, w in edges:
        A[l, i, j] = (A[l, i, j] + w) % p
    return WA(alpha, A, beta)


def _add(p, rows, piv, v):
    """extend the reduced echelon basis (rows, piv) by v; returns the new basis and whether v was independent"""
    no_overflow(p, len(piv))
    r = (v - v[piv] @ rows) % p
    nz = np.flatnonzero(r)
    if len(nz) == 0:
        return rows, piv, False
    c = nz[0]
    r = r * pow(int(r[c]), p - 2, p) % p
    return np.vstack([(rows - np.outer(rows[:, c], r)) % p, r]), piv + [c], True


def _basis(p, starts, step, k, beta=None):
    """BFS over words from the vectors `starts` = [(v, word)]: a reduced echelon basis (rows, pivots) of span{v A_u},
    where step(v, l) = v A_l and k is the number of letters.
    With beta given, stop at the first (shortest) word with v A_u beta != 0 and return it.
    Children of a vector are computed only when popped, to keep the queue small."""
    rows, piv, todo, parent = np.zeros((0, len(starts[0][0])), np.int64), [], deque(starts), {}
    while todo:
        v, u = todo.popleft()
        if v is None:
            v = step(parent[u[:-1]], u[-1])
        if beta is not None and int(v @ beta % p):
            return rows, piv, u
        rows, piv, new = _add(p, rows, piv, v)              # v itself (not its reduction) keeps the word u
        if new:
            parent[u] = v
            todo.extend((None, u + (l,)) for l in range(k))
    return rows, piv, None


def _mul(p, A):
    no_overflow(p, A.shape[1])
    return lambda v, l: v @ A[l] % p


def minimize(p, W):
    """backward then forward reduction; the result has the minimal dimension (Hankel rank mod p).
    Backward goes first: the guessed inputs of fixed qubits at the end of the register make the forward space of a
    loop large (on one benchmark loop 560 states against 45)."""
    W = WA(W.beta, W.A.transpose(0, 2, 1), W.alpha)
    for _ in range(2):
        Fm, piv, _ = _basis(p, [(W.alpha % p, ())], _mul(p, W.A), len(W.A))
        no_overflow(p, Fm.shape[1])
        # coordinates of a vector of span(F) in the reduced basis are its entries at the pivots
        W = WA(W.alpha[piv], np.stack([Fm @ W.A[l] % p for l in range(len(W.A))])[:, :, piv], Fm @ W.beta % p)
        W = WA(W.beta, W.A.transpose(0, 2, 1), W.alpha)          # reverse and repeat
    return WA(W.beta, W.A.transpose(0, 2, 1), W.alpha)


def product(p, W1, W2, terms):
    """automaton with letter l computing sum over (l1, l2) in terms[l] of A1(l1) (x) A2(l2), built on the span reachable
    from alpha1 (x) alpha2 (forward reduction) without forming the Kronecker matrices: a vector is a D1 x D2 matrix V
    and V (A1 (x) A2) = A1^T V A2. Both factors read the same word, so the reachable span is often much smaller than
    D1 D2 (pair of a CX ladder with a loop from qubit 20: 100 of 46^2 = 2116)."""
    D1, D2 = len(W1.alpha), len(W2.alpha)
    no_overflow(p, max(D1, D2))
    step = lambda v, l: sum(W1.A[l1].T @ v.reshape(D1, D2) % p @ W2.A[l2] % p for l1, l2 in terms[l]).ravel() % p
    alpha = np.outer(W1.alpha, W2.alpha).ravel() % p
    Fm, piv, _ = _basis(p, [(alpha, ())], step, len(terms))
    A = np.array([[step(r, l)[piv] for r in Fm] for l in range(len(terms))], np.int64).reshape(len(terms), len(piv), len(piv))
    return WA(alpha[piv], A, (Fm.reshape(-1, D1, D2) @ W2.beta % p) @ W1.beta % p)


def compose(p, W2, W1):
    """automaton of U2 U1 (W1 is applied first): sum over the middle index z of W2(a,x,z) (x) W1(a,z,y)"""
    return product(p, W2, W1, [[(4 * a + 2 * x + z, 4 * a + 2 * z + y) for z in (0, 1)]
                               for a in (0, 1) for x in (0, 1) for y in (0, 1)])


def difference(p, W1, W2):
    """automaton of f1 - f2"""
    D1, D2 = len(W1.alpha), len(W2.alpha)
    A = np.zeros((len(W1.A), D1 + D2, D1 + D2), np.int64)
    A[:, :D1, :D1], A[:, D1:, D1:] = W1.A, W2.A
    return WA(np.concatenate([W1.alpha, W2.alpha]), A, np.concatenate([W1.beta, -W2.beta % p]))


def pair(p, W, Wc):
    """automaton of f(a,x,y) conj f(a,x',y') over letters 16a + 8x + 4y + 2x' + y'; Wc is W built with conj roots"""
    return product(p, W, Wc, [[(4 * a + 2 * x + y, 4 * a + 2 * x2 + y2)]
                              for a in (0, 1) for x in (0, 1) for y in (0, 1) for x2 in (0, 1) for y2 in (0, 1)])


def letters_subset(W, keep):
    return WA(W.alpha, W.A[list(keep)], W.beta)


def nonzero_word(p, W, L=0):
    """shortest word of length >= L with f(w) != 0 mod p, or None. The search starts from a basis of
    {alpha A_u : |u| = L}, built level by level (each level has at most dim W independent vectors)."""
    level = [(W.alpha % p, ())]
    no_overflow(p, len(W.alpha))
    for _ in range(L):
        rows, piv, nxt = np.zeros((0, len(W.alpha)), np.int64), [], []
        for v, u in level:
            for l in range(len(W.A)):
                vl = v @ W.A[l] % p
                rows, piv, new = _add(p, rows, piv, vl)
                if new:
                    nxt.append((vl, u + (l,)))
        level = nxt
    return _basis(p, level, _mul(p, W.A), len(W.A), W.beta)[2] if level else None


def value(p, W, word):
    v = W.alpha
    no_overflow(p, len(W.alpha))
    for l in word:
        v = v @ W.A[l] % p
    return int(v @ W.beta % p)
