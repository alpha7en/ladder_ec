"""Second relations (Table 4, last column): the same pairs with the relation that differs from the task relation
only by the global phase: E2 (up to phase) -> E1 (exact); E5phase -> E5strict (exact on input |0...0>);
E3 (with ancillas) -> E3phase (the same up to phase). Ids, tasks and edits do not change.
The ground truth is in data/benchmark/truth.jsonl (ground_truth.py computes all relations; the truth key of E5strict
is "E5"). The field expect is not set on the new pairs. Output: data/benchmark/pairs_second.jsonl.
Usage: python -I experiments/benchmark/make_second_relations.py"""
import json
import os

BENCH = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "data", "benchmark")
OTHER = {"E2": "E1", "E5phase": "E5strict", "E3": "E3phase"}


def main():
    with open(os.path.join(BENCH, "pairs_second.jsonl"), "w") as f:
        for l in open(os.path.join(BENCH, "pairs.jsonl")):
            p = json.loads(l)
            p.update(label=OTHER[p["label"]], label_task=p["label"])
            p.pop("expect")
            f.write(json.dumps(p, ensure_ascii=False) + "\n")


if __name__ == "__main__":
    main()
