"""Table 4: answers of the method against the ground truth, for the task relation and the second relation
(2200 pairs), and the refusals by kind of Table S2 (Sec. 7.2).
Usage: python -I experiments/tables/table04_answers.py"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import Benchmark, C, KIND, answer, kind, eps_exponent      # noqa: E402

ROWS = [(False, "neq", '"NOT equal" with counterexample'), (False, "refused_sol", "refusal: mutation leaves the class"),
        (False, "timeout", "timeout (300 s)"), (True, "eq_exact", '"equal", exact plan'),
        (True, "eq_eps", '"equal", error <= 2^-40'), (True, "noans", '"no answer", error bound >= 1'),
        (True, "refused_ref", "refusal on the reference"), (True, "refused_sol", "refusal on the solution")]


def count(B, run, label_of):
    cnt, wrong, codes = C(), 0, C()
    for i, p in B.pairs.items():
        L = label_of(i)
        t, o = B.truth_of(i, L), run[i]
        k = kind(answer(o, L))
        if k == "refused":
            side, st = ("ref", o["ref_status"]) if o["ref_status"] != "ladder" else ("sol", o["sol_status"])
            k = "refused_" + side
            code = st.split(":", 1)[1]
            codes[(t, side, KIND.get(code, "?"), code)] += 1
        cnt[(t, k)] += 1
        wrong += (k.startswith("eq") and t is False) or (k == "neq" and t is True)
    return cnt, wrong, codes


def main():
    B = Benchmark()
    main_, second = B.run("method_main.jsonl"), B.run("method_second.jsonl")
    a = count(B, main_, lambda i: B.pairs[i]["label"])
    b = count(B, second, lambda i: B.second[i]["label"])
    print("%-10s %-40s %8s %8s" % ("Truth", "Answer", "Task", "Second"))
    for t, k, name in ROWS:
        print("%-10s %-40s %8d %8d" % ("equal" if t else "not equal", name, a[0][(t, k)], b[0][(t, k)]))
    rest = (set(a[0]) | set(b[0])) - {(t, k) for t, k, _ in ROWS}
    for key in sorted(rest, key=str):
        print("%-10s %-40s %8d %8d" % (key[0], key[1], a[0][key], b[0][key]))
    print("%-51s %8d %8d" % ("wrong answers", a[1], b[1]))
    ex = [eps_exponent(answer(main_[i], B.pairs[i]["label"])) for i in B.equal()
          if kind(answer(main_[i], B.pairs[i]["label"])) == "eq_eps"]
    print("error bounds of the task relation: 2^%d ... 2^%d" % (min(ex), max(ex)))
    print("refusals by kind (task relation): truth, side, kind -> pairs")
    tot = C()
    for (t, side, k, code), v in a[2].items():
        tot[(t, side, k)] += v
    for key in sorted(tot, key=str):
        print("  ", key, tot[key], {c: v for (t, s, kk, c), v in sorted(a[2].items()) if (t, s, kk) == key})


if __name__ == "__main__":
    main()
