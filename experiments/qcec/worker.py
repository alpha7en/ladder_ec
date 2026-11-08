"""Worker process of the QCEC enumeration (experiments/qcec/driver.py). Reads jobs (JSON, one per line) from stdin
and prints one JSON line per job to stdout. It keeps the imported generators between jobs, so that the interpreter
start and the qiskit import are not paid per instance; only building the circuits and checking are timed.

Job: {"ref": {"file", "function"}, "sol": {...}, "n": int, "a": int | null, "crit": "E1|E2|E3|E3phase|E5|E5phase",
      "timeout": float}

QCEC 3.11.0 configuration:
- The ZX checker is off in all modes. With the default portfolio ZX answers "probably_not_equivalent" on equal
  circuits with rotations, and QCEC then reports "no_information" within a fraction of a second although the
  alternating checker proves equality by itself. The alternating and simulation checkers remain, default tolerance.
- E1, E2: qcec.verify(run_zx_checker=False), timeout.
  E1 = answer "equivalent"; E2 = "equivalent" or "equivalent_up_to_global_phase".
- E5, E5phase: all qubits marked as ancillas in |0> (the column U|0...0> is compared), without
  check_partial_equivalence: that mode compares only outcome probabilities and misses phase errors.
  E5 = "equivalent"; E5phase = "equivalent" or "equivalent_up_to_global_phase".
- E3 (exact), E3phase (up to global phase): if the qubit counts differ, the smaller circuit is padded; qubits with
  index >= d, d the smaller qubit count (data first, ancillas last), are marked as ancillas in |0> (not garbage:
  they must return to |0>). If the qubit counts are equal, E3 reduces to E1 (E3phase to E2), which is noted. If the
  counts differ for E1/E2, the check runs as E3 (from E1) or E3phase (from E2), noted in crit_effective.
An exception in the generator itself (undefined at these n, a) is returned as build_error, not as an answer.
Circuits are transpiled to BASIS without optimization (which keeps the global phase), with
qubits_initially_zero=False; the global phase of the circuit is moved into gates (explicit_phase), otherwise QCEC
does not see it and cannot tell E1 from E2.
"""
import importlib.util
import inspect
import json
import os
import resource
import sys
import time
import traceback

from mqt import qcec
from mqt.core import load
from qiskit import QuantumCircuit, transpile

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
BASIS = ["h", "x", "y", "z", "s", "sdg", "t", "tdg", "rx", "ry", "rz", "p", "cx", "cz", "cp", "ccx", "swap"]
_mods = {}


def get_fn(ref):
    path = ref["file"]
    if not os.path.isabs(path):
        path = os.path.join(ROOT, path)
    if path not in _mods:
        d = os.path.dirname(path)
        if d not in sys.path:
            sys.path.insert(0, d)
        spec = importlib.util.spec_from_file_location("m%d" % len(_mods), path)
        m = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(m)
        _mods[path] = m
    return getattr(_mods[path], ref["function"])


def build(ref, n, a):
    fn = get_fn(ref)
    params = list(inspect.signature(fn).parameters)
    args = [n] if a is None or len(params) < 2 else [n, a]
    qc = fn(*args)
    # qubits_initially_zero=False is required: by default Qiskit 2.x treats untouched qubits as |0> and uses
    # them as clean ancillas when decomposing mcx; the circuit is then correct only on part of the inputs
    # (false "not equal" answers otherwise).
    tq = transpile(qc, basis_gates=BASIS, optimization_level=0, qubits_initially_zero=False)
    return explicit_phase(tq)


def explicit_phase(qc):
    """QCEC 3.11 ignores the global_phase of a circuit: it calls Z and e^{i pi/2} RZ(pi) (exactly equal)
    "equivalent_up_to_global_phase". For E1 and E5 (exact) to make sense, the phase is moved into gates:
    P(t) X P(t) X = e^{it} I on qubit 0."""
    t = float(qc.global_phase)
    if abs(t) < 1e-15:
        return qc
    out = qc.copy()
    out.global_phase = 0
    out.p(t, 0); out.x(0); out.p(t, 0); out.x(0)
    return out


def pad(qc, m):
    if qc.num_qubits >= m:
        return qc
    out = QuantumCircuit(m, global_phase=qc.global_phase)
    out.compose(qc, range(qc.num_qubits), inplace=True)
    return out


def run(job):
    n, a, crit = job["n"], job.get("a"), job["crit"]
    t0 = time.perf_counter()
    try:
        c1, c2 = build(job["ref"], n, a), build(job["sol"], n, a)
    except Exception as e:
        return {"build_error": repr(e)[:200], "build_s": round(time.perf_counter() - t0, 4)}
    note = None
    opts = {"timeout": float(job["timeout"]), "run_zx_checker": False}
    eff = crit
    if crit in ("E1", "E2") and c1.num_qubits != c2.num_qubits:
        eff = "E3" if crit == "E1" else "E3phase"
        note = f"qubit counts {c1.num_qubits} and {c2.num_qubits}: checked as E3"
    if eff in ("E3", "E3phase"):
        d = min(c1.num_qubits, c2.num_qubits)
        if job.get("e3_free") is not None:
            # restricted input of the task (only qubit 0 free, the others |0>)
            d = int(job["e3_free"])
            note = f"E3 as a restricted input: qubits 0..{d - 1} free, the others ancillas in |0>"
        m = max(c1.num_qubits, c2.num_qubits)
        c1, c2 = pad(c1, m), pad(c2, m)
        q1, q2 = load(c1), load(c2)
        if m > d:
            for q in (q1, q2):
                q.set_circuit_qubits_ancillary(d, m - 1)
        elif note is None:
            note = "E3 with equal qubit counts: checked as E1 or E2"
    else:
        if c1.num_qubits != c2.num_qubits:
            return {"error": f"different qubit counts {c1.num_qubits} and {c2.num_qubits} for {crit}",
                    "build_s": round(time.perf_counter() - t0, 4)}
        q1, q2 = load(c1), load(c2)
        if crit in ("E5", "E5phase"):
            for q in (q1, q2):
                q.set_circuit_qubits_ancillary(0, q.num_qubits - 1)
            # check_partial_equivalence is NOT used: in QCEC it is equality of measurement probability
            # distributions, not of state vectors (a GHZ mutant with an extra S would pass).
    t1 = time.perf_counter()
    res = qcec.verify(q1, q2, **opts)
    t2 = time.perf_counter()
    eq = str(res.equivalence).split(".")[-1]
    exact = eq == "equivalent"
    upto = eq in ("equivalent", "equivalent_up_to_global_phase")
    if eq == "no_information" or eq == "probably_equivalent" or eq == "probably_not_equivalent":
        verdict = "unknown" if eq == "no_information" else eq
    elif eff in ("E1", "E5", "E3"):
        verdict = "equal" if exact else ("not_equal" if eq in ("not_equivalent", "equivalent_up_to_global_phase", "equivalent_up_to_phase") else eq)
    else:
        verdict = "equal" if upto else ("not_equal" if eq in ("not_equivalent", "equivalent_up_to_phase") else eq)
    timed_out = (t2 - t1) >= float(job["timeout"]) * 0.98 and verdict not in ("equal", "not_equal")
    out = {"qcec_result": eq, "crit_effective": eff,
           "checkers": [[c.get("checker"), c.get("equivalence")] for c in res.json().get("checkers", [])], "verdict": "timeout" if timed_out else verdict,
           "qubits": c1.num_qubits, "gates_ref": c1.size(), "gates_sol": c2.size(),
           "build_s": round(t1 - t0, 4), "check_s": round(t2 - t1, 4),
           "maxrss_mb": round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 1)}
    if note:
        out["note"] = note
    if out["verdict"] == "not_equal" and c1.num_qubits <= 12:
        out["qiskit_check"] = qiskit_check(c1, c2, eff)
    return out


def qiskit_check(c1, c2, crit):
    """Recheck of "not equal" with the Qiskit matrix or state vector (float64), tolerance 1e-8: tells a real
    counterexample from a false answer of QCEC."""
    import numpy as np
    from qiskit.quantum_info import Operator, Statevector
    if crit in ("E5", "E5phase"):
        s1, s2 = Statevector(c1).data, Statevector(c2).data
        if crit == "E5":
            return "not_equal" if np.abs(s1 - s2).max() > 1e-8 else "equal"
        return "not_equal" if abs(abs(np.vdot(s1, s2)) - 1) > 1e-8 else "equal"
    U1, U2 = Operator(c1).data, Operator(c2).data
    if crit in ("E3", "E3phase"):
        return "skipped (E3)"
    if crit == "E1":
        return "not_equal" if np.abs(U1 - U2).max() > 1e-8 else "equal"
    k = np.unravel_index(np.argmax(np.abs(U1)), U1.shape)
    c = U2[k] / U1[k]
    return "not_equal" if abs(abs(c) - 1) > 1e-8 or np.abs(U2 - c * U1).max() > 1e-8 else "equal"


def main():
    print(json.dumps({"ready": True, "qcec": qcec.__version__ if hasattr(qcec, "__version__") else "3.11.0"}), flush=True)
    for line in sys.stdin:
        if not line.strip():
            continue
        job = json.loads(line)
        try:
            out = run(job)
        except Exception as e:
            out = {"error": repr(e)[:300], "trace": traceback.format_exc()[-600:]}
        print(json.dumps(out, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
