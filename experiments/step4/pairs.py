"""Pairs of the step-4 experiment: (name, reference, solution, relation, K); the threshold of the pair is n_0 = K + 2."""
PAIRS = [(f"{fam}_k{K}", f"{fam}_ref", f"{fam}_k{K}", "E2", K) for fam, Ks in
         (("cx", (5, 8, 11, 14, 20)), ("rzz", (5, 8, 11, 14, 20)), ("bv", (5, 8, 11))) for K in Ks]
