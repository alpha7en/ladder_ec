"""Check of the dimension bound of Theorem 1 (Sec. 5.2, Table 1) on every distinct loop of the benchmark ladders
(Sec. 7.2, Online Resource 1 S4): the number of states of the constructed automaton, modulo the first prime
p = 1 (mod m), must not exceed the bound. Also counts the loops where the bound is attained ("tight").
Usage: python -I experiments/checks/theorem1_bound.py"""
import collections
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from ladder_shape import loops                         # noqa: E402
from ladder_ec.ladders import segment_wa               # noqa: E402
from ladder_ec.field import Field, field_order         # noqa: E402
from ladder_ec.plans import primes_in_order            # noqa: E402
from ladder_ec.decide import angles                    # noqa: E402


def main():
    stats, fails = collections.Counter(), []
    for fn, file, fam, sg, sh in loops():
        m = field_order(angles(fam))
        M = len(segment_wa(Field(next(primes_in_order(m)), m), sg, fam.extra).alpha)
        stats["loops"] += 1
        stats["tight"] += M == sh["bound"]
        if M <= sh["bound"]:
            stats["bound_ok"] += 1
        else:
            stats["bound_FAIL"] += 1
            fails.append((fn, M, sh))
    print(dict(stats))
    for f in fails:
        print("FAIL", *f)
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
