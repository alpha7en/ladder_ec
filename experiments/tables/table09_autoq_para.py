"""Table 9: the method and AutoQ-Para on shared families (experiments/autoq_para/): median of three runs, s.
AutoQ-Para time: wall time of /usr/bin/time; for the GHZ row the experiment comparing both programs on |0>.
Usage: python -I experiments/tables/table09_autoq_para.py"""
import json, os, re, statistics, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import jl, run_path      # noqa: E402

AUTOQ_OF = {"ghz_ladder_vs_elevator_E5": "ladder_vs_elevator_on_zero"}    # the experiment timed for this row


def autoq_times(path):
    t = {}
    for line in open(path):
        m = re.search(r'Running experiment "([^"]+)".*Status : (\w+).*wall=([\d.]+)s', line)
        if m:
            t.setdefault(m.group(1), []).append((m.group(2), float(m.group(3))))
    return t


def main():
    rows = jl(run_path("autoq_para/ours.jsonl"))
    truth = json.load(open(run_path("autoq_para/truth.json")))
    at = autoq_times(run_path("autoq_para/autoq_runs.txt"))
    print("%-28s %-10s %-34s %10s %10s" % ("Pair", "Relation", "Answer of the method", "AutoQ-Para", "Method"))
    for r in rows:
        exp = AUTOQ_OF.get(r["pair"], r["autoq"])
        a = "%.2f (%s)" % (statistics.median(w for _, w in at[exp]), at[exp][0][0]) if exp in at else "-"
        v = r["verdict"] if isinstance(r["verdict"], str) else r["verdict"][r["relation"]]
        ok = (v.startswith("equal") == truth[r["pair"]]["truth"]) if not v.startswith("refused") else None
        print("%-28s %-10s %-34s %10s %10.2f  correct=%s" % (r["pair"], r["relation"], v, a, r["median_s"], ok))


if __name__ == "__main__":
    main()
