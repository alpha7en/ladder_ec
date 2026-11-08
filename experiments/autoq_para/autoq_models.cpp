// Models of ladder families for AutoQ-Para (Table 9): the GHZ ladder and related experiments.
// Built together with the AutoQ-Para sources; experiments.cpp is included for the functions of its Table 1.
#define main autoq_experiments_main
#include "experiments.cpp"
#undef main

#include <cstring>
#include <functional>

using DLF = std::vector<Def_Linear_Form>;
using ACN = Algebraic_Complex_Number;

static const Internal_Symbol SYM_X  = 0;  // ordinary qubit
static const Internal_Symbol SYM_XP = 1;  // second-to-last qubit (terminal symbol of the ladder)

// |0...0> on n >= 2 qubits; the second-to-last qubit is marked by the symbol x'.
SWTA make_zero_state() {
    SWTA::Metadata md = {.number_of_internal_symbols = 2, .number_of_colors = 1};
    SWTA::Transition_Builder b(md);
    State q = 0, q_last = 1, leaf = 2, bot = 3;
    Def_Coef one(ACN::ONE()), zero(ACN::ZERO());
    auto add = [&](State s, Internal_Symbol sym, const Def_Linear_Form& l, const Def_Linear_Form& r) {
        DLF L{l}, R{r};
        b.add_transition(s, sym, 0, synthetize_swta_transition(L, R));
    };
    add(q, SYM_X, one * q, zero * bot);
    add(q, SYM_XP, one * q_last, zero * bot);
    add(q_last, SYM_X, one * leaf, zero * bot);
    b.add_bot_state_transitions(bot);
    Bit_Set leaves(4, {leaf, bot});
    return SWTA(b.build(4), {q}, leaves);
}

// All basis states (as HAMILTONIAN_ALL_BASIS, but through our helper).
SWTA make_all_basis() { return get_predefined_swta(Predefined_SWTA_Names::HAMILTONIAN_ALL_BASIS); }

// (|0...0> + |1...1>)/sqrt(2) on n >= 2 qubits.
SWTA make_ghz_state() {
    SWTA::Metadata md = {.number_of_internal_symbols = 2, .number_of_colors = 1};
    SWTA::Transition_Builder b(md);
    State r = 0, z = 1, o = 2, z_last = 3, o_last = 4, leaf = 5, bot = 6;
    Def_Coef one(ACN::ONE()), zero(ACN::ZERO()), c(ACN::ONE_OVER_SQRT2());
    auto add = [&](State s, Internal_Symbol sym, const Def_Linear_Form& l, const Def_Linear_Form& rr) {
        DLF L{l}, R{rr};
        b.add_transition(s, sym, 0, synthetize_swta_transition(L, R));
    };
    add(r, SYM_X, c * z, c * o);
    add(r, SYM_XP, c * z_last, c * o_last);
    add(z, SYM_X, one * z, zero * bot);
    add(z, SYM_XP, one * z_last, zero * bot);
    add(o, SYM_X, zero * bot, one * o);
    add(o, SYM_XP, zero * bot, one * o_last);
    add(z_last, SYM_X, one * leaf, zero * bot);
    add(o_last, SYM_X, zero * bot, one * leaf);
    b.add_bot_state_transitions(bot);
    Bit_Set leaves(7, {leaf, bot});
    return SWTA(b.build(7), {r}, leaves);
}

// Single-qubit gate [[a,b],[c,d]] on the first qubit, identity on the others.
WTT make_first_qubit_gate(const ACN& a, const ACN& bb, const ACN& cc, const ACN& d) {
    SWTA::Metadata md = {.number_of_internal_symbols = 2, .number_of_colors = 1};
    WTT_Builder b(md);
    State s0 = 0, id = 1;
    b.mark_state_initial(s0);
    b.mark_state_final(id);
    Def_Coef A(a), B(bb), C(cc), D(d), one(ACN::ONE());
    for (Internal_Symbol sym : {SYM_X, SYM_XP}) {
        DLF L{A * id * Subtree_Tag::LEFT}, R{D * id * Subtree_Tag::RIGHT};
        if (!bb.is_zero()) L.push_back(B * id * Subtree_Tag::RIGHT);
        if (!cc.is_zero()) R.push_back(C * id * Subtree_Tag::LEFT);
        b.add_transition(s0, sym, synthetize_wtt_transition(L, R));
        DLF Li{one * id * Subtree_Tag::LEFT}, Ri{one * id * Subtree_Tag::RIGHT};
        b.add_transition(id, sym, synthetize_wtt_transition(Li, Ri));
    }
    return b.build(2);
}

WTT make_h_first() {
    return make_first_qubit_gate(ACN::ONE_OVER_SQRT2(), ACN::ONE_OVER_SQRT2(), ACN::ONE_OVER_SQRT2(), -ACN::ONE_OVER_SQRT2());
}

// Two-qubit box: a single-qubit gate on one of the two qubits (pos=0 first, 1 second).
WTT make_box_1q(int pos, const ACN& a, const ACN& bb, const ACN& cc, const ACN& d) {
    SWTA::Metadata md = {.number_of_internal_symbols = 1, .number_of_colors = 1};
    WTT_Builder b(md);
    State s0 = 0, s1 = 1, leaf = 2;
    b.mark_state_initial(s0);
    b.mark_state_final(leaf);
    Def_Coef A(a), B(bb), C(cc), D(d), one(ACN::ONE());
    auto gate = [&](State from, State to) {
        DLF L{A * to * Subtree_Tag::LEFT}, R{D * to * Subtree_Tag::RIGHT};
        if (!bb.is_zero()) L.push_back(B * to * Subtree_Tag::RIGHT);
        if (!cc.is_zero()) R.push_back(C * to * Subtree_Tag::LEFT);
        b.add_transition(from, 0, synthetize_wtt_transition(L, R));
    };
    auto idt = [&](State from, State to) {
        DLF L{one * to * Subtree_Tag::LEFT}, R{one * to * Subtree_Tag::RIGHT};
        b.add_transition(from, 0, synthetize_wtt_transition(L, R));
    };
    if (pos == 0) { gate(s0, s1); idt(s1, leaf); } else { idt(s0, s1); gate(s1, leaf); }
    return b.build(3);
}

// CX with control on the second qubit of the box and target on the first (reversed direction).
WTT make_box_cx_reversed() {
    SWTA::Metadata md = {.number_of_internal_symbols = 1, .number_of_colors = 1};
    WTT_Builder b(md);
    // Output (y0,y1) = (x0 xor x1, x1). Read y0: the left subtree (y0=0) takes x0=y1, the right one (y0=1) takes x0=1-y1.
    // States: s0 -> (keep, flip); keep: y1 from x(L/R) as usual; flip: left and right for x0 crossed.
    // Done with second-level linear forms: s0 reads x0, but the choice of the x0 subtree depends on y1,
    // so it is encoded as a sum of two terms (the zeroL/zeroR construction of Sec. 6.1 of the AutoQ-Para paper).
    State s0 = 0, pL = 1, pR = 2, leaf = 3;
    b.mark_state_initial(s0);
    b.mark_state_final(leaf);
    Def_Coef one(ACN::ONE());
    // y0 = 0: (x0=0,x1=0)->(0,0), (x0=1,x1=1)->(0,1). Left output of s0 = pL(L-branch of x0) "+" pR(R-branch of x0).
    {
        DLF L{one * pL * Subtree_Tag::LEFT, one * pR * Subtree_Tag::RIGHT};
        DLF R{one * pR * Subtree_Tag::LEFT, one * pL * Subtree_Tag::RIGHT};
        b.add_transition(s0, 0, synthetize_wtt_transition(L, R));
    }
    // pL lets through only y1=0 from x1=0, pR only y1=1 from x1=1.
    Def_Coef zero(ACN::ZERO());
    {
        DLF L{one * leaf * Subtree_Tag::LEFT}, R{zero * leaf * Subtree_Tag::RIGHT};
        b.add_transition(pL, 0, synthetize_wtt_transition(L, R));
    }
    {
        DLF L{zero * leaf * Subtree_Tag::LEFT}, R{one * leaf * Subtree_Tag::RIGHT};
        b.add_transition(pR, 0, synthetize_wtt_transition(L, R));
    }
    return b.build(4);
}

WTT ladder_stage(WTT box) { return perform_staircase_construction(box, {0, 0}, 1, SYM_XP); }

WTT cx_box() { return get_predefined_wtt(Predefined_WTT_Names::HAMILTONIAN_BC_CNOT, {}); }

// CX as H.CZ.H on the target.
WTT cx_box_via_cz() {
    ACN s = ACN::ONE_OVER_SQRT2();
    WTT h2 = make_box_1q(1, s, s, s, -s);
    // CZ: phase -1 on |11>.
    SWTA::Metadata md = {.number_of_internal_symbols = 1, .number_of_colors = 1};
    WTT_Builder b(md);
    State s0 = 0, i1 = 1, z1 = 2, leaf = 3;
    b.mark_state_initial(s0);
    b.mark_state_final(leaf);
    Def_Coef one(ACN::ONE()), m1(-ACN::ONE());
    auto tr = [&](State from, Def_Coef cl, State tl, Def_Coef cr, State trr) {
        DLF L{cl * tl * Subtree_Tag::LEFT}, R{cr * trr * Subtree_Tag::RIGHT};
        b.add_transition(from, 0, synthetize_wtt_transition(L, R));
    };
    tr(s0, one, i1, one, z1);
    tr(i1, one, leaf, one, leaf);
    tr(z1, one, leaf, m1, leaf);
    WTT cz = b.build(4);
    WTT t1 = compose_wtts_sequentially(h2, cz);
    return compose_wtts_sequentially(t1, h2);
}

// Hand-written transducer of the elevator: CX(0, i) for all i >= 1 (y_0 = x_0, y_i = x_i xor x_0).
// Algorithm 1 of AutoQ-Para does not build it (the control qubit is fixed, the box does not shift); written by hand.
WTT make_elevator() {
    SWTA::Metadata md = {.number_of_internal_symbols = 2, .number_of_colors = 1};
    WTT_Builder b(md);
    State s0 = 0, c0 = 1, c1 = 2;
    b.mark_state_initial(s0);
    b.mark_state_final(c0);
    b.mark_state_final(c1);
    Def_Coef one(ACN::ONE());
    for (Internal_Symbol sym : {SYM_X, SYM_XP}) {
        { DLF L{one * c0 * Subtree_Tag::LEFT}, R{one * c1 * Subtree_Tag::RIGHT}; b.add_transition(s0, sym, synthetize_wtt_transition(L, R)); }
        { DLF L{one * c0 * Subtree_Tag::LEFT}, R{one * c0 * Subtree_Tag::RIGHT}; b.add_transition(c0, sym, synthetize_wtt_transition(L, R)); }
        { DLF L{one * c1 * Subtree_Tag::RIGHT}, R{one * c1 * Subtree_Tag::LEFT}; b.add_transition(c1, sym, synthetize_wtt_transition(L, R)); }
    }
    return b.build(3);
}

Verification_Status run_program_vs_swta(const SWTA& pre, std::vector<WTT> stages, const SWTA& post) {
    SWTA cur = pre;
    for (auto& w : stages) cur = apply_wtt_to_swta(cur, w);
    return are_two_swtas_color_equivalent(cur, post) ? SUCCESS : FAILURE;
}

Verification_Status run_two_programs(const SWTA& pre, std::vector<WTT> p1, std::vector<WTT> p2) {
    SWTA a = pre, b = pre;
    for (auto& w : p1) a = apply_wtt_to_swta(a, w);
    for (auto& w : p2) b = apply_wtt_to_swta(b, w);
    return are_two_swtas_color_equivalent(a, b) ? SUCCESS : FAILURE;
}

int main(int argc, char** argv) {
    std::map<std::string, std::function<Verification_Status()>> exps = {
        {"bv", verify_bv},
        {"grover", verify_grover},
        {"adder", verify_adder},
        {"ecc", verify_ecc},
        {"hamiltonian", verify_hamiltonian_simulation},
        // E5: H(0), then the ladder CX(i,i+1) from |0..0> gives GHZ.
        {"ghz_ladder_e5", [] { return run_program_vs_swta(make_zero_state(), {make_h_first(), ladder_stage(cx_box())}, make_ghz_state()); }},
        // Mutant: ladder with reversed CX (control on i+1) from |0..0> does NOT give GHZ. Expected FAIL.
        {"ghz_ladder_rev_e5_mutant", [] { return run_program_vs_swta(make_zero_state(), {make_h_first(), ladder_stage(make_box_cx_reversed())}, make_ghz_state()); }},
        // Mutant: without H. Expected FAIL.
        {"ghz_ladder_noH_e5_mutant", [] { return run_program_vs_swta(make_zero_state(), {ladder_stage(cx_box())}, make_ghz_state()); }},
        // E1 on all basis states: CX ladder against the H.CZ.H ladder. Expected SUCCESS.
        {"ladder_cx_vs_czh_allbasis", [] { return run_two_programs(make_all_basis(), {ladder_stage(cx_box())}, {ladder_stage(cx_box_via_cz())}); }},
        // E1 on all basis states: CX ladder against the ladder of reversed CX. Expected FAIL.
        {"ladder_cx_vs_rev_allbasis_mutant", [] { return run_two_programs(make_all_basis(), {ladder_stage(cx_box())}, {ladder_stage(make_box_cx_reversed())}); }},
        // Model check: reversed CX = (H x H) CX (H x H) inside the box. Expected SUCCESS.
        {"ladder_rev_vs_hcxh_allbasis", [] {
            ACN s = ACN::ONE_OVER_SQRT2();
            WTT h0 = make_box_1q(0, s, s, s, -s), h1 = make_box_1q(1, s, s, s, -s);
            WTT hh = compose_wtts_sequentially(h0, h1);
            WTT cx = cx_box();
            WTT t = compose_wtts_sequentially(hh, cx);
            WTT rev2 = compose_wtts_sequentially(t, hh);
            return run_two_programs(make_all_basis(), {ladder_stage(make_box_cx_reversed())}, {ladder_stage(rev2)}); }},
        // E5: H(0) and the elevator from |0..0> give GHZ. Expected SUCCESS.
        {"ghz_elevator_e5", [] { return run_program_vs_swta(make_zero_state(), {make_h_first(), make_elevator()}, make_ghz_state()); }},
        // E5 by comparing two programs: ladder against elevator on |0..0>. Expected SUCCESS.
        {"ladder_vs_elevator_on_zero", [] { return run_two_programs(make_zero_state(), {make_h_first(), ladder_stage(cx_box())}, {make_h_first(), make_elevator()}); }},
        // E1 on all basis states: ladder against elevator. Unitarily different, expected FAIL.
        {"ladder_vs_elevator_allbasis", [] { return run_two_programs(make_all_basis(), {ladder_stage(cx_box())}, {make_elevator()}); }},
    };
    if (argc < 2) {
        for (auto& [k, v] : exps) std::cout << k << "\n";
        return 0;
    }
    std::string name = argv[1];
    if (!exps.contains(name)) { std::cerr << "unknown\n"; return 2; }
    RUN_EXPERIMENT(name, exps[name]);
    return 0;
}
