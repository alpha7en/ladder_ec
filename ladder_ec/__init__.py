"""Parameterized equivalence checking of quantum circuit generators via weighted automata.

check(reference, solution, relation) decides a relation between two Python generator functions
solve(n[, a]) -> QuantumCircuit for all n and all bit strings a (Algorithm 1 of the paper).
"""
from .check import check, decide_core, RELATIONS

__version__ = "1.0.0"
__all__ = ["check", "decide_core", "RELATIONS", "__version__"]
