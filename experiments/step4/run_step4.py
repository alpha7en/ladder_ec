"""Step-4 experiment (Table 8): step 4 of the method (sizes n < n_0, dense matrices modulo p, instances.check with
exact=True) against the QCEC enumeration of the same sizes and all a (experiments/qcec/driver.py, mode A, bound
n_0 - 1; 60 s per instance, 600 s per pair).
  ours  decide_core on every pair with a limit of 600 s: automata time, step-4 time, answer
        -> data/runs/step4_ours.jsonl
  qcec  QCEC enumeration at n = 1..n_0 - 1 -> data/runs/step4_qcec.jsonl (input step4_qcec_input.jsonl)
Usage: python -I experiments/step4/run_step4.py ours|qcec"""
import json, os, signal, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path[:0] = [HERE, ROOT, os.path.join(ROOT, "experiments", "qcec")]
from pairs import PAIRS                       # noqa: E402
RUNS = os.path.join(ROOT, "data", "runs")
G = os.path.join(HERE, "gens.py")
REL = "experiments/step4/gens.py"


class Timeout(Exception):
    pass


def ours():
    from ladder_ec.check import decide_core
    from ladder_ec.translator import analyze
    signal.signal(signal.SIGALRM, lambda *_: (_ for _ in ()).throw(Timeout()))
    decide_core((G, "cx_ref"), (G, "cx_ref"), ["E2"])          # warm-up
    with open(os.path.join(RUNS, "step4_ours.jsonl"), "w") as f:
        for name, r, s, rel, K in PAIRS:
            n0 = max(analyze(G, r)["n_min"], analyze(G, s)["n_min"])
            t = time.time()
            signal.alarm(600)
            try:
                res = decide_core((G, r), (G, s), [rel])
            except Timeout:
                res = {"verdicts": "timeout 600 s"}
            except MemoryError:
                res = {"verdicts": "MemoryError"}
            finally:
                signal.alarm(0)
            row = dict(pair=name, K=K, n0=n0, verdict=res.get("verdicts"), time_s=round(time.time() - t, 3),
                       automata_s=res.get("automata_s"), quick_search_s=res.get("quick_search_s"),
                       step4_s=res.get("step4_s"), step4=res.get("step4"))
            f.write(json.dumps(row, ensure_ascii=False, default=str) + "\n")
            f.flush()
            print(row, flush=True)


def qcec():
    import driver as D
    from ladder_ec.translator import analyze
    D.PAIR_LIMIT = 600.0
    D.CONFIG["limit_pair_s"] = 600.0
    D.MEM_LIMIT = 3 << 30
    D.CONFIG["tool"] += "; run from experiments/step4/run_step4.py (sizes n < n_0, pair limit 600 s)"
    inp = os.path.join(RUNS, "step4_qcec_input.jsonl")
    with open(inp, "w") as f:
        for name, r, s, rel, K in PAIRS:
            n0 = max(analyze(G, r)["n_min"], analyze(G, s)["n_min"])
            f.write(json.dumps(dict(id=name, task=name, label=rel, ref=dict(file=REL, function=r),
                                    sol=dict(file=REL, function=s), origin="review", cutoff={rel: n0 - 1},
                                    n_min=1)) + "\n")
    sys.argv = [sys.argv[0], "--mode", "A", "--pairs", inp, "--jobs", "1", "--out", os.path.join(RUNS, "step4_qcec.jsonl"),
                "--inst", os.path.join(RUNS, "step4_qcec_instances.jsonl")]
    D.main()


if __name__ == "__main__":
    {"ours": ours, "qcec": qcec}[sys.argv[1]]()
