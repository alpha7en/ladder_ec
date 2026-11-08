"""Command line: python -m ladder_ec REFERENCE.py:function SOLUTION.py:function [--relation R] [--seed S] [--json]"""
import argparse
import json
import sys

from .check import check, RELATIONS


def main(argv=None):
    ap = argparse.ArgumentParser(prog="ladder_ec", description="Equivalence of two quantum circuit generators "
                                 "for all n and all bit strings a.")
    ap.add_argument("reference", help="path.py:function of the reference generator")
    ap.add_argument("solution", help="path.py:function of the solution generator")
    ap.add_argument("--relation", default="eq_phase", choices=sorted(RELATIONS),
                    help="eq (=), eq_phase (= up to phase), eq0 / eq0_phase (on |0>), "
                         "eq_anc / eq_anc_phase (with returned ancillas); default eq_phase")
    ap.add_argument("--seed", type=int, help="seed of the probabilistic plan (written to the output)")
    ap.add_argument("--json", action="store_true", help="print the full result as JSON")
    a = ap.parse_args(argv)
    r = check(a.reference, a.solution, a.relation, seed=a.seed)
    if a.json:
        print(json.dumps(r, indent=1, default=str))
    else:
        print(r["answer"])
        if r["counterexample"]:
            print("counterexample at n = %d (found by %s)" % (r["counterexample"]["n"], r["counterexample"]["found_by"]))
        d = r["details"]
        if d["ref_status"] != "ladder" or d["sol_status"] != "ladder":
            print("reference: %s; solution: %s" % (d["ref_status"], d["sol_status"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
