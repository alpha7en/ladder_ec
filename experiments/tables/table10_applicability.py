"""Table 10: generators outside the benchmark (RQ4, Sec. 7.5): how many functions the translator accepts and why it
refuses. "With n": by the signature. The kind of a refusal is the column "Kind" of Table S2 (common.KIND):
out of class, out of fragment, or kind undetermined (both).
MQT Bench and QuanBench refuse everything; for them the class was decided by reading the code
(data/external/mqtbench/classification.csv; QuanBench: the 14 functions with a size parameter call the library
Fourier transform). Usage: python -I experiments/tables/table10_applicability.py"""
import csv, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import C, KIND, ROOT, jl, run_path      # noqa: E402


def kinds(codes):
    c = C(KIND.get(x, "?" + x) for x in codes)
    return c["class"], c["fragment"], c["both"], {k: v for k, v in c.items() if k.startswith("?")}


def line(name, total, with_n, accepted, cls, frag, both, extra=""):
    print("%-38s %9s %7s %9s %9s %9s %9s %s" % (name, total, with_n, accepted, cls, frag, both, extra))


def main():
    line("Source", "Functions", "With n", "Accepted", "Out of class", "Out of fragment", "Undetermined")
    mq = jl(run_path("mqtbench_strip_defaults.jsonl"))
    cl = list(csv.DictReader(open(os.path.join(ROOT, "data", "external", "mqtbench", "classification.csv"))))
    assert len(mq) == len(cl)
    line("MQT Bench", len(mq), len(mq), sum(r["status"] == "ladder" for r in mq),
         sum(r["class"] == "out of class" for r in cl), sum(r["class"] != "out of class" for r in cl), 0)
    own = list(csv.DictReader(open(run_path("own_generators_refusals.csv"))))
    k = kinds(r["reason"] for r in own if r["status"] != "ladder")
    line("own generators for other tasks", len(own), len(own), sum(r["status"] == "ladder" for r in own), *k[:3], k[3] or "")
    sel = jl(run_path("katas_selected.jsonl"))
    kt = [r for r in jl(run_path("katas_translator.jsonl")) if r["group"] == "family"]
    k = kinds(r["B"]["reason"] for r in kt if r["B"]["status"] != "ladder")
    line("QuantumKatas, references", len(sel), len(kt), sum(r["B"]["status"] == "ladder" for r in kt), *k[:3], k[3] or "")
    llm = jl(run_path("llm_katas.jsonl"))
    k = kinds(r["B"] for r in llm if r["B"] != "ladder")
    line("QuantumKatas, language-model solutions", len(llm), len(llm), sum(r["B"] == "ladder" for r in llm), *k[:3],
         k[3] or "")
    qb = [r for r in jl(run_path("quanbench_translate.jsonl")) if not r["source"].startswith("QuanBench+")]
    wn = [r for r in qb if r["size_param"]]
    line("QuanBench", len(qb), len(wn), sum(r["translator"] == "ladder" for r in qb), len(wn), 0, 0)
    print("translator codes on QuanBench:", dict(C(r["translator"] for r in qb)))


if __name__ == "__main__":
    main()
