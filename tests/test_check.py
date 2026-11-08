"""End-to-end: the GHZ example of the README and the public entry function."""
import os

import pytest

from conftest import GEN
from ladder_ec import check

G = os.path.join(GEN, "ghz.py")


def test_ghz_same_state():
    r = check(G + ":ghz_ref", G + ":ghz_star", "eq0", seed=1)
    assert r["answer"].startswith("equal") and r["counterexample"] is None


def test_ghz_different_unitary():
    r = check(G + ":ghz_ref", G + ":ghz_star", "eq", seed=1)
    assert r["answer"] == "NOT equal" and r["counterexample"]["n"] == 3


def test_ghz_phase_flip():
    assert check(G + ":ghz_ref", G + ":ghz_phase", "eq0", seed=1)["answer"] == "NOT equal"
    assert check(G + ":ghz_ref", G + ":ghz_phase", "eq0_phase", seed=1)["answer"] == "NOT equal"


def test_unknown_relation():
    with pytest.raises(ValueError):
        check(G + ":ghz_ref", G + ":ghz_ref", "equivalent")
