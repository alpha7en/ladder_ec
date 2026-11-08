"""Table 8: constructed correct solutions with threshold n_0 = K + 2 (experiments/step4/): time of step 4 of the
method (sizes 1 <= n < n_0, exact, Proposition 2) and of the QCEC enumeration of the same sizes, s. "---": the
run did not finish within 600 s and 8 GB. Usage: python -I experiments/tables/table08_step4.py"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import jl, run_path      # noqa: E402

FAMILIES = [("cx", "CX chain"), ("rzz", "RZZ chain"), ("bv", "BV oracle")]


def main():
    ours = {r["pair"]: r for r in jl(run_path("step4_ours.jsonl"))}
    qcec = {r["id"]: r for r in jl(run_path("step4_qcec.jsonl"))}
    n0s = sorted({r["n0"] for r in ours.values()})
    print("%-12s %-7s" % ("Reference", "") + "".join("%9s" % ("n0=%d" % n) for n in n0s))
    for fam, name in FAMILIES:
        by_n0 = {r["n0"]: r for p, r in ours.items() if p.startswith(fam + "_")}
        cell = lambda r: "---" if r["step4_s"] is None else "%.3g" % r["step4_s"]
        print("%-12s %-7s" % (name, "method") + "".join("%9s" % (cell(by_n0[n]) if n in by_n0 else "")
                                                         for n in n0s))
        qc = lambda p: "%.3g" % qcec[p]["wall_s"] if qcec[p]["verdict"] == "equal" else "---"
        print("%-12s %-7s" % ("", "QCEC") + "".join("%9s" % (qc(by_n0[n]["pair"]) if n in by_n0 else "")
                                                     for n in n0s))
    print("answers of the method:", {p: r["verdict"] for p, r in ours.items()})


if __name__ == "__main__":
    main()
