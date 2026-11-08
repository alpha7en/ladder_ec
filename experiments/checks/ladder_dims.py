"""Data of Fig. S2: the bound of Theorem 1 against the number of states of the constructed automaton and against
the smallest dimension, for every distinct loop of the benchmark ladders (Sec. 7.2, Online Resource 1 S4).
The smallest dimension is the dimension after minimization modulo the first three primes p = 1 (mod m); the largest
of the three is taken, since the rank modulo p is at most the rank over the field.
Output: data/runs/ladder_dims.csv and counts on screen.
Usage: python -I experiments/checks/ladder_dims.py"""
import csv
import itertools
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from ladder_shape import ROOT, loops                   # noqa: E402
from ladder_ec.ladders import segment_wa               # noqa: E402
from ladder_ec.field import Field, field_order         # noqa: E402
from ladder_ec.plans import primes_in_order            # noqa: E402
from ladder_ec.decide import angles                    # noqa: E402
from ladder_ec.automata import minimize                # noqa: E402


def main():
    rows = []
    for fn, file, fam, sg, sh in loops():
        m = field_order(angles(fam))
        Ms, rs = [], []
        for p in itertools.islice(primes_in_order(m), 3):
            W = segment_wa(Field(p, m), sg, fam.extra)
            Ms.append(len(W.alpha))
            rs.append(len(minimize(p, W).alpha))
        rows.append(dict(function=fn, file=file, s=sh["s"], w=sh["w"], S=sh["S"], E=sh["E"], h=sh["h"], b=sh["b"],
                         i0=sh["i0"], K0=sh["K0"], m=m, bound=sh["bound"], main_term=sh["main_term"], D0=sh["D0"],
                         states=max(Ms), min_dim=max(rs), min_dim_per_prime="/".join(map(str, rs)),
                         states_per_prime="/".join(map(str, Ms))))
    with open(os.path.join(ROOT, "data", "runs", "ladder_dims.csv"), "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=list(rows[0]))
        wr.writeheader()
        wr.writerows(rows)
    print("ladders", len(rows))
    print("bound >= states:", sum(r["bound"] >= r["states"] for r in rows),
          " bound == states:", sum(r["bound"] == r["states"] for r in rows))
    print("states >= min_dim:", sum(r["states"] >= r["min_dim"] for r in rows),
          " states == min_dim:", sum(r["states"] == r["min_dim"] for r in rows))
    print("bound == min_dim:", sum(r["bound"] == r["min_dim"] for r in rows))
    print("min_dim differs between primes:", sum(len(set(r["min_dim_per_prime"].split("/"))) > 1 for r in rows))


if __name__ == "__main__":
    main()
