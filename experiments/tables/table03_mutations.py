"""Table 3: mutation kinds and their pair counts (benchmark without no-op mutations, Sec. 7.1).
Usage: python -I experiments/tables/table03_mutations.py"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import Benchmark, C

# mutation kind codes of the benchmark generator -> rows of Table 3
ROW = {"G": "gate", "F": "phase", "C": "loop", "I": "index", "Bb": "condition on bits a", "R": "branch on n",
       "U": "angle", "A": "ancilla"}
ORDER = ["gate", "loop", "index", "condition on bits a", "branch on n", "angle", "ancilla", "phase"]


def row_of(kind):
    if kind.startswith("style"):
        return None
    head = kind.split("_")[0].rstrip("0123456789")
    return ROW[head]


def main():
    B = Benchmark()
    c = C(row_of(p["kind"]) for p in B.pairs.values())
    print("pairs %d, tasks %d; relations %s" % (len(B.pairs), len({p["task"] for p in B.pairs.values()}),
                                                 dict(C(p["label"] for p in B.pairs.values()))))
    print("correct solutions: %d reference-reference, %d rewritings; mutants %d" % (
        sum(p["kind"] == "style:self" for p in B.pairs.values()),
        sum(p["kind"].startswith("style") and p["kind"] != "style:self" for p in B.pairs.values()),
        sum(v for k, v in c.items() if k)))
    print("%-22s %5s" % ("Kind", "Pairs"))
    for r in ORDER:
        print("%-22s %5d" % (r, c[r]))


if __name__ == "__main__":
    main()
