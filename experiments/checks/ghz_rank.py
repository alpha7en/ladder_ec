"""The GHZ chain H(0); for i := 1 to n-1 do CX(i-1, i): automaton of Theorem 1, its minimization and the Hankel rank
(the example with window w = 2 of Sec. 5).

Usage: python -I experiments/checks/ghz_rank.py

What is checked.
1. The automaton of Theorem 1 (ladder_ec/ladders.py: fixed_wa, loop_wa, family_wa) modulo two primes p = 1 (mod 8):
   the number of states of the transition graphs of the segments and the dimension after minimization (F3).
2. Definition 5 (Sec. 5.1): the value of the automaton on every word of length N = n, 2 <= n <= 7, equals <x|U(n)|y>
   modulo p; U(n) is computed separately, in floating point, from the circuit (qubit 0 is the first tensor factor).
3. Rank of blocks of the Hankel matrix H[u, v] = f(uv) over words u, v with |uv| in [2, 7] (n <= 7, values defined by
   the circuit): exactly over Q (the entries sqrt2 * <x|U|y> are integers 0, +-1), by Gaussian elimination in fractions.
4. The explicit automaton with three states q_i, q_0, q_1 (the figure): agreement with U(n) on all words of length 1..7.
Letters (a, x, y) with index 4a + 2x + y, as in ladder_ec; the bit a has no effect (the program has no conditions).
"""
import sys, itertools
from fractions import Fraction
import os
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import numpy as np
from ladder_ec.field import Field
from ladder_ec.ladders import Gate, Fixed, Loop, Family, c, n, i, segment_wa, family_wa
from ladder_ec.automata import minimize, value

NMAX = 7
fam = Family([Fixed([Gate("h", [c(0)])]),
              Loop(c(1), n(0), 1, [Gate("cx", [i(-1), i(0)])])], 0, 2)   # range(1, n): i = 1..n-1; n_min = 2


# ---------- the circuit in floating point ----------
S2 = np.sqrt(2.0)
Hm = np.array([[1, 1], [1, -1]]) / S2


def U_of(nq):
    N = nq
    U = np.eye(2 ** N)
    ops = [np.eye(2)] * N
    ops[0] = Hm
    M = ops[0]
    for o in ops[1:]:
        M = np.kron(M, o)
    U = M @ U
    for ii in range(1, N):                 # CX(i-1, i), i = 1..n-1, in the order of application
        P = np.zeros((2 ** N, 2 ** N))
        for k in range(2 ** N):
            b = [(k >> (N - 1 - j)) & 1 for j in range(N)]
            if b[ii - 1]:
                b[ii] ^= 1
            P[sum(v << (N - 1 - j) for j, v in enumerate(b)), k] = 1
        U = P @ U
    return U


Us = {nq: U_of(nq) for nq in range(1, NMAX + 1)}


def f_exact(word):
    """sqrt2 * <x|U|y> as an integer; word is a tuple of letters (x, y)"""
    N = len(word)
    xi = int("".join(str(x) for x, _ in word), 2)
    yi = int("".join(str(y) for _, y in word), 2)
    v = Us[N][xi, yi] * S2
    r = round(v)
    assert abs(v - r) < 1e-9
    return r


def formula(word):                          # closed form: <x0|H|y0> prod [x_j = x_{j-1} xor y_j]
    (x0, y0) = word[0]
    v = (-1) ** (x0 * y0)                    # sqrt2 <x0|H|y0>
    for j in range(1, len(word)):
        if word[j][0] != word[j - 1][0] ^ word[j][1]:
            return 0
    return v


def rank_Q(rows):
    M = [[Fraction(v) for v in r] for r in rows]
    rk, cols = 0, len(M[0]) if M else 0
    for col in range(cols):
        piv = next((r for r in range(rk, len(M)) if M[r][col] != 0), None)
        if piv is None:
            continue
        M[rk], M[piv] = M[piv], M[rk]
        for r in range(len(M)):
            if r != rk and M[r][col] != 0:
                q = M[r][col] / M[rk][col]
                M[r] = [a - q * b for a, b in zip(M[r], M[rk])]
        rk += 1
    return rk


XY = [(x, y) for x in (0, 1) for y in (0, 1)]


def words(L):
    return list(itertools.product(XY, repeat=L))


print("GHZ chain: H(0); for i := 1 to n-1 step 1 do CX(i-1, i) od;  e = 0, n_min = 2")

# 2. closed form against the circuit
mx = max(abs(f_exact(wd) - formula(wd)) for L in range(1, NMAX + 1) for wd in words(L))
print("closed form <x0|H|y0> prod [x_j = x_{j-1} xor y_j] against the circuit, n = 1..%d: max difference (x sqrt2) = %d" % (NMAX, mx))

# 1. automaton of Theorem 1 and minimization
for p in (16777289, 16777337):              # the first primes p = 1 (mod 8), p >= 2^24 (as the exact plan, Sec. 6.2)
    F = Field(p, 8)
    segs = [segment_wa(F, seg, fam.extra) for seg in fam.segments]
    raw = [len(W.alpha) for W in segs]
    mins = [len(minimize(p, W).alpha) for W in segs]
    W = family_wa(F, fam)
    # Definition 5: the value of the automaton on words of length n = 2..7 equals <x|U|y> mod p
    s2inv = pow((pow(F.r, 1, p) + pow(F.r, 7, p)) % p, p - 2, p)      # 1/sqrt2 = 1/(zeta_8 + zeta_8^{-1})
    bad = 0
    total = 0
    for L in range(2, NMAX + 1):
        for wd in words(L):
            for a in (0, 1):                 # the bit a has no effect: both letters are checked
                lw = [4 * a + 2 * x + y for x, y in wd]
                total += 1
                if value(p, W, lw) != f_exact(wd) * s2inv % p:
                    bad += 1
    print("p = %d: states of the transition graph of H(0): %d, of the ladder: %d; after minimization: %d and %d; "
          "automaton of the program after minimization: %d; mismatches with <x|U|y> on %d words of length 2..%d: %d"
          % (p, raw[0], raw[1], mins[0], mins[1], len(W.alpha), total, NMAX, bad))

# 3. ranks of Hankel blocks (exactly over Q)
def block(ulens, vlens):
    rows = [wu for k in ulens for wu in words(k)]
    cols = [wv for k in vlens for wv in words(k)]
    return [[f_exact(wu + wv) for wv in cols] for wu in rows], len(rows), len(cols)

for ulens, vlens in [((1, 2, 3), (1, 2, 3, 4)), ((0, 1, 2, 3), (2, 3, 4)), ((0, 1, 2), (2, 3, 4, 5)),
                     ((0, 1, 2, 3), (2, 3, 4))]:
    assert min(ulens) + min(vlens) >= 2 and max(ulens) + max(vlens) <= NMAX
    M, r, cc = block(ulens, vlens)
    print("Hankel block |u| in %s, |v| in %s (%d x %d): rank over Q = %d" % (ulens, vlens, r, cc, rank_Q(M)))

# the column U|0^n> (letters with y = 0): for comparison with "two states"
rows = [wu for k in (0, 1, 2, 3) for wu in itertools.product([(0, 0), (1, 0)], repeat=k)]
cols = [wv for k in (2, 3) for wv in itertools.product([(0, 0), (1, 0)], repeat=k)]
print("series of the column U|0^n> (y = 0), |u| <= 3, |v| in (2, 3): rank over Q = %d" % rank_Q([[f_exact(u + v) for v in cols] for u in rows]))

# 4. explicit automaton of the figure: states q_i, q_0, q_1; weights multiplied by sqrt2 only on transitions from q_i
def explicit(wd):
    st = {"i": 1}
    for x, y in wd:
        nxt = {}
        for q, wgt in st.items():
            if q == "i":
                nxt[x] = nxt.get(x, 0) + wgt * (-1) ** (x * y)        # weight sqrt2 <x|H|y>
            elif x == q ^ y:
                nxt[x] = nxt.get(x, 0) + wgt
        st = nxt
    return st.get(0, 0) + st.get(1, 0)                                  # final weights: q_0, q_1 -> 1, q_i -> 0

bad = sum(explicit(wd) != f_exact(wd) for L in range(1, NMAX + 1) for wd in words(L))
print("explicit automaton (q_i, q_0, q_1) against the circuit on all words of length 1..%d: %d mismatches" % (NMAX, bad))

# 5. why the automaton f_P of Theorem 1 has 4 states after minimization and not 3: on words of length < n_min + e = 2
# it need not agree with U (the ladder is empty at n = 1, Lemma S.3), and its series is a different series on short words.
p = 16777289
F = Field(p, 8)
W = family_wa(F, fam)
s2inv = pow((F.r + pow(F.r, 7, p)) % p, p - 2, p)
short = {}
for L in (0, 1):
    for wd in words(L):
        lw = [2 * x + y for x, y in wd]
        v = value(p, W, lw)
        short[wd] = 0 if v == 0 else (1 if v == s2inv else (-1 if v == (-s2inv) % p else None))
print("automaton f_P on words of length 0 and 1 (x sqrt2):", {("".join("%d%d" % b for b in k) or "eps"): v for k, v in short.items()})

def fP(wd):
    return short[wd] if len(wd) < 2 else f_exact(wd)

for ulens, vlens in [((0, 1, 2, 3), (0, 1, 2, 3, 4)), ((0, 1, 2), (0, 1, 2, 3, 4, 5))]:
    rows = [wu for k in ulens for wu in words(k)]
    cols = [wv for k in vlens for wv in words(k)]
    print("series f_P: Hankel block |u| in %s, |v| in %s: rank over Q = %d" % (ulens, vlens, rank_Q([[fP(u + v) for v in cols] for u in rows])))
