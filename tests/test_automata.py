"""Automaton of a ladder program against the matrix of the expanded circuit (Thm. 1): every entry
<y|U(n, a)|x> equals the value of the automaton on the word of (a_j, x_j, y_j)."""
import itertools
import random
from fractions import Fraction

import pytest

import conftest  # noqa: F401
from ladder_ec.field import Field, random_prime, field_order
from ladder_ec.ladders import Family, Fixed, Loop, Gate, c, n, i, expand, matrix, family_wa
from ladder_ec.automata import value

FAMILIES = {
    "ghz_ladder": Family([Fixed([Gate("h", [c(0)])]), Loop(c(0), n(-1), 1, [Gate("cx", [i(0), i(1)])])], 0, 2),
    "ghz_star": Family([Fixed([Gate("h", [c(0)])]), Loop(c(1), n(0), 1, [Gate("cx", [c(0), i(0)])])], 0, 2),
    "ghz_down": Family([Fixed([Gate("h", [n(-1)])]), Loop(n(-1), c(0), -1, [Gate("cx", [i(0), i(-1)])])], 0, 2),
    "bv_anc": Family([Fixed([Gate("x", [n(0)]), Gate("h", [n(0)])]), Loop(c(0), n(0), 1, [Gate("h", [i(0)])]),
                      Loop(c(0), n(0), 1, [Gate("cx", [i(0), n(0)], cond=((i(0), 1),))]),
                      Loop(c(0), n(0), 1, [Gate("h", [i(0)])])], 1, 1),
    "bv_phase": Family([Loop(c(0), n(0), 1, [Gate("z", [i(0)], cond=((i(0), 1),))])], 0, 1),
    "win3": Family([Loop(c(2), n(-2), 1, [Gate("ccx", [i(0), i(1), i(-1)]), Gate("t", [i(1)]),
                                         Gate("cx", [c(0), i(0)], cond=((i(1), 0),)),
                                         Gate("rz", [n(-1)], Fraction(3, 8))])], 0, 5),
    "down2": Family([Fixed([Gate("h", [c(0)]), Gate("cz", [c(1), n(-1)])]),
                     Loop(n(-2), c(1), -1, [Gate("ch", [i(0), i(1)]), Gate("cx", [c(0), i(0)], cond=((i(1), 1),))])],
                    0, 4),
    "brick": Family([Fixed([Gate("h", [c(0)]), Gate("h", [c(1)])]), Loop(c(0), n(-1), 2, [Gate("cx", [i(0), i(1)])]),
                     Loop(c(1), n(-1), 2, [Gate("cz", [i(0), i(1)]), Gate("t", [i(1)], cond=((i(1), 1),))])], 0, 3),
    "gap3": Family([Loop(c(0), n(-1), 3, [Gate("h", [i(0)]), Gate("ch", [i(0), i(1)]), Gate("t", [i(1)])])], 0, 2),
    "down2s": Family([Fixed([Gate("h", [c(0)]), Gate("x", [n(-1)])]),
                      Loop(n(-2), c(1), -2, [Gate("cx", [c(0), i(0)]), Gate("swap", [i(0), i(-1)]),
                                             Gate("cz", [n(-1), i(-1)], cond=((i(0), 0),))])], 0, 4),
}


@pytest.mark.parametrize("name", sorted(FAMILIES))
def test_automaton_matches_matrix(name):
    fam = FAMILIES[name]
    rng = random.Random(3)
    angles = [g.angle for s in fam.segments for g in (s.gates if isinstance(s, Fixed) else s.body)
              if g.angle is not None]
    F = Field(random_prime(field_order(angles), rng), field_order(angles))
    W = family_wa(F, fam)
    for nn in range(fam.n_min, fam.n_min + 3):
        N = nn + fam.extra
        for _ in range(3):
            a = [rng.randint(0, 1) for _ in range(N)]
            M = matrix(F, expand(fam, nn, a), N)
            for xy in itertools.product((0, 1), repeat=2 * N):
                x, y = xy[:N], xy[N:]
                assert int(M[xy]) == value(F.p, W, [4 * a[j] + 2 * x[j] + y[j] for j in range(N)])
