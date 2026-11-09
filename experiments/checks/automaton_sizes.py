"""Size of the automaton of Theorem 1 before and after minimization for the examples of Sec. 5.1: the reference R of
Example 1 and the GHZ chain of Example 2.

Usage: python -I experiments/checks/automaton_sizes.py

Computed modulo the two first primes p = 1 (mod 8), p >= 2^24:
1. the number of states of the automaton of each segment by Theorem 1 (transition graph of fixed_wa / loop_wa,
   function explore) and the dimension after minimization (F3);
2. the dimension of the product by Lemma 1 without minimization: the product of the numbers of states of the segments;
3. the dimension of the automaton of the program after minimization (family_wa: minimization of the segments and of
   each product);
4. the bound of Theorem 1 for the ladder: 4^h 2^b ((4^w - 1)/3 + max(0, s - w)) + D_0,
   D_0 = i_0 8^|S| 2^|E| + 4^|E| (s (K_0 + 1) + s (s - 1)/2) (Online Resource 1, S1.3); the parameters w, h, S, E, i_0,
   K_0 are those computed by loop_wa; b is the number of bits of a the automaton stores (here 0: the bit a_i of the
   condition is read in the same letter where the body is applied);
5. check: the value of the automaton of the program on words of length N = n + e, n_min <= n <= 4, equals
   <x|U(n,a)|y> modulo p (U from ladders.expand and ladders.matrix, the bit a_j is taken from the letter of qubit j).
Programs:
  R (Example 1): e = 1; X(n); H(n); for i := 0 to n-1 step 1 do if a[i] = 1 then CX(i, n) fi od; H(n); X(n).
  GHZ (Example 2): e = 0; H(0); for i := 1 to n-1 step 1 do CX(i-1, i) od.
"""
import sys, itertools
import os
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from ladder_ec.field import Field
from ladder_ec.ladders import Gate, Fixed, Loop, Family, c, n, i, segment_wa, family_wa, expand, matrix
from ladder_ec.automata import minimize, value

PRIMES = (16777289, 16777337)          # the first primes p = 1 (mod 8), p >= 2^24 (as in fig_ghz_rank.py)

R = Family([Fixed([Gate("x", [n(0)]), Gate("h", [n(0)])]),
            Loop(c(0), n(0), 1, [Gate("cx", [i(0), n(0)], None, ((i(0), 1),))]),     # range(0, n): i = 0..n-1
            Fixed([Gate("h", [n(0)]), Gate("x", [n(0)])])], 1, 1)
GHZ = Family([Fixed([Gate("h", [c(0)])]),
              Loop(c(1), n(0), 1, [Gate("cx", [i(-1), i(0)])])], 0, 2)               # range(1, n): i = 1..n-1


def ladder_params(seg, e):
    """w, h, |S|, |E|, s, i_0, K_0 as in ladders.loop_wa (positive step)"""
    refs = [q for g in seg.body for q in g.qubits] + [r for g in seg.body for r, _ in g.cond]
    offs = [k for kind, k in refs if kind == "i"]
    lo, w = min(offs), max(offs) - min(offs) + 1
    i0, dd = seg.start[1] + lo, seg.stop[1] + lo - e
    K0 = 1 - dd - w
    S = sorted({k for kind, k in refs if kind == "c"})
    E = sorted({k for kind, k in refs if kind == "n"})
    return dict(w=w, S=len(S), E=len(E), h=len(S) + len(E), s=seg.step, i0=i0, K0=K0)


def bound(P, b=0):
    w, h, S, E, s, i0, K0 = (P[k] for k in ("w", "h", "S", "E", "s", "i0", "K0"))
    first = 4 ** h * 2 ** b * ((4 ** w - 1) // 3 + max(0, s - w))
    D0 = i0 * 8 ** S * 2 ** E + 4 ** E * (s * (K0 + 1) + s * (s - 1) // 2)
    return first, D0


def check(F, fam, W, nmax=4):
    """mismatches between the value of the automaton and <x|U(n,a)|y> on all words of length n + e, n_min <= n <= nmax"""
    p, e = F.p, fam.extra
    bad = total = 0
    for nn in range(fam.n_min, nmax + 1):
        N = nn + e
        for a in itertools.product((0, 1), repeat=N):
            M = matrix(F, expand(fam, nn, list(a)), N)
            for x in itertools.product((0, 1), repeat=N):
                for y in itertools.product((0, 1), repeat=N):
                    word = [4 * a[j] + 2 * x[j] + y[j] for j in range(N)]
                    total += 1
                    bad += value(p, W, word) != int(M[x + y]) % p
    return bad, total


for name, fam in (("R (Example 1)", R), ("GHZ (Example 2)", GHZ)):
    print("== %s: e = %d, n_min = %d" % (name, fam.extra, fam.n_min))
    for seg in fam.segments:
        if isinstance(seg, Loop):
            P = ladder_params(seg, fam.extra)
            first, D0 = bound(P)
            print("   ladder: w = %(w)d, h = %(h)d (|S| = %(S)d, |E| = %(E)d), s = %(s)d, i_0 = %(i0)d, K_0 = %(K0)d, b = 0" % P)
            print("   bound of Theorem 1 for the ladder: %d + D_0 = %d + %d = %d" % (first, first, D0, first + D0))
    for p in PRIMES:
        F = Field(p, 8)
        segs = [segment_wa(F, seg, fam.extra) for seg in fam.segments]
        raw = [len(W.alpha) for W in segs]
        mins = [len(minimize(p, W).alpha) for W in segs]
        prod_raw = 1
        for r in raw:
            prod_raw *= r
        prod_min = 1
        for r in mins:
            prod_min *= r
        W = family_wa(F, fam)
        bad, total = check(F, fam, W)
        kinds = ["ladder" if isinstance(s, Loop) else "fixed" for s in fam.segments]
        print("   p = %d: segments (%s): states %s, after minimization %s;" % (p, ", ".join(kinds), raw, mins))
        print("      product by Lemma 1 without minimization: %s = %d; product of the minimized segments: %d;"
              % (" x ".join(map(str, raw)), prod_raw, prod_min))
        print("      automaton of the program after minimization (family_wa): %d; mismatches with <x|U|y> on %d words: %d"
              % (len(W.alpha), total, bad))


# Why R has 4 states after minimization and 3 in Fig. 5a. The automaton A_R of the figure (q_e, q_o, q_f) defines the
# family (Definition 5: agreement on words of length N = n + e >= n_min + e = 2). The automaton of Theorem 1 has to
# agree with U only on these words; on words of length 0 and 1 (n < n_min) its series f_P is different, and
# minimization (F3) keeps the series f_P on all words, so it gives 4. The same holds for GHZ (fig_ghz_rank.py): 3 for
# the smallest automaton of the family, 4 for f_P.
def A_R(word):
    """explicit automaton of Fig. 5a: letter (a, x, y); x != y gives the zero matrix; q_e -> q_f and q_o -> q_f on any letter"""
    st = {"e": 1}
    for a, x, y in word:
        if x != y:
            return 0
        nxt = {}
        for q, wgt in st.items():
            if q == "f":
                continue
            par = (q == "o") ^ (a * x)
            nxt["o" if par else "e"] = nxt.get("o" if par else "e", 0) + wgt
            nxt["f"] = nxt.get("f", 0) + wgt * (1 if q == "e" else (-1) ** (1 - x))
        st = nxt
    return st.get("f", 0)


p = PRIMES[0]
F = Field(p, 8)
W = family_wa(F, R)
bad = total = 0
for nn in range(1, 5):
    N = nn + 1
    for a in itertools.product((0, 1), repeat=N):
        M = matrix(F, expand(R, nn, list(a)), N)
        for x in itertools.product((0, 1), repeat=N):
            for y in itertools.product((0, 1), repeat=N):
                total += 1
                bad += A_R(list(zip(a, x, y))) % p != int(M[x + y]) % p
print("== R: explicit automaton A_R with three states (Fig. 5a) against <x|U(n,a)|y>, n = 1..4: %d mismatches of %d" % (bad, total))
diff = []
for L in (0, 1):
    for wd in itertools.product(range(8), repeat=L):
        fa = A_R([((l >> 2), (l >> 1) & 1, l & 1) for l in wd]) % p
        fp = value(p, W, list(wd))
        if fa != fp:
            diff.append(("".join(map(str, wd)) or "eps", fp if fp < p // 2 else fp - p, fa if fa < p // 2 else fa - p))
print("   words of length 0 and 1 where the series f_P (automaton of Theorem 1) and A_R differ (word, f_P, A_R):", diff)
