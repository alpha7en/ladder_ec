"""Pairs of the AutoQ-Para comparison (Table 9): (name, reference, solution, relation, expected answer, AutoQ-Para
experiment, what AutoQ-Para checks). E1: matrix equality; E5strict: state equality from |0...0>."""
PAIRS = [
    ("ghz_ladder_vs_elevator_E5", "ghz_elevator", "ghz_ladder", "E5strict", True, "ghz_ladder_e5; ghz_elevator_e5; ladder_vs_elevator_on_zero", "ladder and elevator from |0> give GHZ"),
    ("ghz_noH_E5", "ghz_ladder", "ghz_ladder_noh", "E5strict", False, "ghz_ladder_noH_e5_mutant", "mutant without H does not give GHZ"),
    ("ghz_rev_E5", "ghz_ladder", "ghz_ladder_rev", "E5strict", False, "ghz_ladder_rev_e5_mutant", "mutant with reversed CX does not give GHZ"),
    ("cx_ladder_vs_czh_E1", "cx_ladder", "cx_ladder_czh", "E1", True, "ladder_cx_vs_czh_allbasis", "all basis inputs"),
    ("cx_ladder_vs_rev_E1", "cx_ladder", "cx_ladder_rev", "E1", False, "ladder_cx_vs_rev_allbasis_mutant", "all basis inputs"),
    ("cx_rev_vs_hcxh_E1", "cx_ladder_rev", "cx_ladder_rev_hcxh", "E1", True, "ladder_rev_vs_hcxh_allbasis", "all basis inputs"),
    ("cx_ladder_vs_elevator_E1", "cx_elevator", "cx_ladder", "E1", False, "ladder_vs_elevator_allbasis", "all basis inputs"),
    ("bv10_E5", "bv10_post", "bv10", "E5strict", True, "bv", "precondition and postcondition, hidden string (10)*"),
    ("bv10_mut_E5", "bv10_post", "bv10_mut", "E5strict", False, None, None),
    ("bv_bits_E5", "bv_bits_post", "bv_bits", "E5strict", True, None, "not in AutoQ-Para: any hidden string"),
    ("adder_uma3_vs_uma2_E1", "adder_uma2", "adder_uma3", "E1", True, "adder", "precondition and postcondition (addition)"),
    ("adder_mut_E1", "adder_uma2", "adder_mut", "E1", False, None, None),
    ("syndrome_two_passes_E1", "syndrome", "syndrome_two_passes", "E1", True, "ecc", "precondition and postcondition (syndrome of one error)"),
    ("syndrome_mut_E1", "syndrome", "syndrome_mut", "E1", False, None, None),
    ("heis_opt_vs_naive_E1", "heis_naive", "heis_opt", "E1", True, "hamiltonian", "equivalence of two circuits"),
    ("heis_mut_E1", "heis_naive", "heis_mut", "E1", False, None, None),
    ("grover_ideal_self_E1", "grover_ideal", "grover_ideal", "E1", True, "grover", "equivalence of two circuits"),
]
