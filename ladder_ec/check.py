"""Algorithm 1 of the paper on two generator functions: translator, steps 1-4, answer (Sec. 6, Table 2).

Order of the stages:
  translator   both functions -> ladder programs, or a refusal with a reason code (Thm. 3, Table S2);
  step 1       different ancilla counts: no unitary relation can hold, any (n, a) is a witness;
  quick search direct simulation modulo a random prime up to 13 qubits (instances.py); a difference found
               there is an exact counterexample and saves the automata;
  steps 2-3    automata and zeroness modulo the primes of a plan (decide.py, plans.py, Thm. 2);
  step 4       the sizes 1 <= n < n_0 on the code itself, exactly (Proposition 2).

Relation keys (decide.py): E1, E2, E5strict, E5phase, E3, E3phase; the public names of RELATIONS map to them.
"""
import math
import secrets
import time
from fractions import Fraction

from . import instances as SN
from .translator import analyze
from .ladders import Gate, Fixed, Loop, Family
from .decide import decide, NO_EPS
from .automata import USED

RELATIONS = {"eq": "E1", "eq_phase": "E2", "eq0": "E5strict", "eq0_phase": "E5phase",
             "eq_anc": "E3", "eq_anc_phase": "E3phase"}
KINDS = {"E1": "E1", "E2": "E2", "E5strict": "E5", "E5phase": "E5p", "E3": "E3", "E3phase": "E3p"}
PRE_QUBITS, PRE_MAX_A = 13, 8               # quick search: direct simulation up to 13 qubits, 8 bit strings per n


def family(d):
    """ladder/v1 description (translator output) -> Family of ladders.py"""
    ref = lambda r: ({"const": "c", "n": "n", "i": "i"}[r["ref"]], r["off"])
    def gate(g):
        cond = tuple((ref(b["index"]), b["value"]) for b in g.get("cond", []))      # nested if: all must hold
        return Gate(g["op"], [ref(q) for q in g["qubits"]], g.get("angle") and Fraction(g["angle"]["pi"]), cond)
    segs = [Fixed([gate(g) for g in s["gates"]]) if s["kind"] == "fixed" else
            Loop(ref(s["range"]["start"]), ref(s["range"]["stop"]), s["range"]["step"], [gate(g) for g in s["body"]])
            for s in d["segments"]]
    return Family(segs, d["num_qubits"]["off"], d["n_min"])


def _status(d):
    return d["status"] if d["status"] == "ladder" else "refused:" + d["reason"]["code"]


def decide_core(ref, sol, kinds, ff=False, seed=None, quick_search=True):
    """Run Algorithm 1 for the relation keys `kinds` on ref, sol = (path, function name).

    ff: restricted input (qubit 0 arbitrary, the other qubits |0>), decide.first_free.
    seed: of the probabilistic plan; drawn afresh when None and written to the output, so that a rerun with it repeats
    the run. quick_search=False is the ablation of Sec. 7 without the quick search.
    Returns a dict with the statuses, the answer per key ("verdicts"), the counterexamples and the timings."""
    seed = secrets.randbits(32) if seed is None else seed
    dR, dG = analyze(*ref), analyze(*sol)
    ours = {"ref_status": _status(dR), "sol_status": _status(dG), "seed": seed}
    if dR["status"] != "ladder" or dG["status"] != "ladder":
        # A solution that reads a bit a[k] with k >= n is outside the domain where the reference is defined
        # (Sec. 3); with a passed as an integer the read gives 0, so both functions are run by the quick search and
        # a difference found there is an exact counterexample.
        run_it = dR["status"] == "ladder" and ours["sol_status"] == "refused:bit_index_out_of_range"
        ours["verdicts"] = "refused (translator)"
        if run_it and quick_search:
            t0 = time.time()
            try:
                pre, _ = SN.check(ref, sol, PRE_QUBITS + 1 - dR["num_qubits"]["off"], kinds, primes=1,
                                  max_a=PRE_MAX_A, free_first=ff)
            except Exception as ex:                  # outside what the simulation can run (angles, parameters)
                pre = {}
                ours["quick_search_error"] = repr(ex)[:200]
            ours["quick_search_s"] = round(time.time() - t0, 3)
            found = {k: v for k, v in pre.items() if v[0] != "equal"}
            if found and len(found) == len(kinds):
                ours["verdicts"] = {k: "NOT equal" if v == "NOT equal" else "reference does not return its ancillas"
                                    for k, (v, i) in found.items()}
                ours["counterexamples"] = {k: {"n": i["n"], "found_by": "quick_search_bit_index"}
                                           for k, (v, i) in found.items() if v == "NOT equal"}
        return ours
    nmin, eR, eG = max(dR["n_min"], dG["n_min"]), dR["num_qubits"]["off"], dG["num_qubits"]["off"]
    # step 1: unitaries on different numbers of qubits are never equal; any (n, a) of the reference is a witness
    out = {k: ("NOT equal", dict(n=dR["n_min"], found_by="ancilla_counts")) for k in kinds
           if not k.startswith("E3") and eR != eG}
    t0 = time.time()
    pre, _ = SN.check(ref, sol, PRE_QUBITS + 1 - max(eR, eG), [k for k in kinds if k not in out], primes=1,
                      max_a=PRE_MAX_A, free_first=ff) if quick_search else ({}, [])
    out.update({k: ("NOT equal" if v == "NOT equal" else "reference does not return its ancillas", dict(info, found_by="quick_search"))
                for k, (v, info) in pre.items() if v != "equal"})
    ours["quick_search_s"] = round(time.time() - t0, 3)
    G, R = family(dG), family(dR)
    primes, cuts = {}, {}
    for k in kinds:                                                  # steps 2-3
        if k in out:
            continue
        v, info = decide(G, R, KINDS[k], seed, free_first=ff)
        info["found_by"] = "automata"
        primes[k] = info.get("primes")
        if "n_cut" in info:
            cuts[k] = info["n_cut"]
        if v == "no answer" and info.get("reason") == NO_EPS:
            v = "no answer (error bound >= 1)"
        elif v == "no answer":
            v = "no answer (no primes)"
        elif v == "undefined":
            v = "reference does not return its ancillas"
        elif v == "equal" and info["error"] == 0:
            v = "equal, exact"
        elif v == "equal":                                           # Table 2: certified only below eps_0 = 2^-40
            e = math.floor(math.log2(info["error"]))
            v = "equal, error <= 2^%d" % (e + 1) if e < -40 else "equal (NOT certified, error bound 2^%d)" % (e + 1)
        out[k] = ("NOT equal" if v == "not equal" else v, info)
    ours["automata_s"] = round(time.time() - t0 - ours["quick_search_s"], 3)
    ours["n0"] = nmin                         # step 4 checks the sizes 1 <= n < n0 on the code
    left = [k for k in kinds if out[k][0].startswith("equal")]
    if left:
        t1 = time.time()
        small, ns = SN.check(ref, sol, nmin, left, free_first=ff, exact=True)
        ours["step4_s"] = round(time.time() - t1, 3)
        ours["step4"] = {"n": ns, "verdicts": {k: v for k, (v, i) in small.items()}}
        for k in left:
            if small[k][0] != "equal":
                out[k] = ("NOT equal" if small[k][0] == "NOT equal" else "reference does not return its ancillas",
                          dict(small[k][1], below_n_min=True, found_by="step4"))
    ours["used"] = dict(USED, step4_max_p=SN.USED_P[0])     # running maxima in this process (Sec. 6.2)
    ours["primes"] = primes
    ours["n_cut"] = cuts                      # enumerating n_0 <= n <= n_cut also decides (Corollary 1)
    ours["verdicts"] = {k: v for k, (v, i) in out.items()}
    ours["counterexamples"] = {k: {"n": i["n"], "found_by": i["found_by"]} for k, (v, i) in out.items()
                               if v == "NOT equal"}
    return ours


def _fn(x):
    """(path, function) or "path:function" -> (path, function)"""
    if isinstance(x, str):
        path, _, name = x.rpartition(":")
        return path, name
    return tuple(x)


def check(reference, solution, relation="eq_phase", seed=None, restricted_input=False):
    """Decide `relation` between two generators for all n and all bit strings a.

    reference, solution: (path, function name) or "path.py:function"; the function is `solve(n)` or
    `solve(n, a)` returning a QuantumCircuit. relation: a key of RELATIONS (eq, eq_phase, eq0, eq0_phase, eq_anc,
    eq_anc_phase) or an internal key (E1, E2, E5strict, E5phase, E3, E3phase).
    Returns a dict: "answer" (Table 2), "counterexample" ({"n", "found_by"} or None) and "details", the full
    output of decide_core."""
    key = RELATIONS.get(relation, relation)
    if key not in KINDS:
        raise ValueError("unknown relation %r; use one of %s" % (relation, ", ".join(RELATIONS)))
    d = decide_core(_fn(reference), _fn(solution), [key], ff=restricted_input, seed=seed)
    v = d["verdicts"]
    return {"answer": v if isinstance(v, str) else v[key],
            "counterexample": d.get("counterexamples", {}).get(key), "details": d}
