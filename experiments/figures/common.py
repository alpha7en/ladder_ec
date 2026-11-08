"""Shared helpers of the figure scripts: data paths, the deterministic horizontal jitter, optional plotting."""
import hashlib
import os

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RUNS = os.path.join(ROOT, "data", "runs")
BENCH = os.path.join(ROOT, "data", "benchmark")
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")


def jitter(key, salt="", width=0.04):
    """u in [-width, width], deterministic in the key: plotted values are spread by 10^u (drawing only)"""
    return (int(hashlib.sha256((key + salt).encode()).hexdigest()[:8], 16) / 0xFFFFFFFF * 2 - 1) * width


def pyplot():
    """matplotlib.pyplot with a non-interactive backend, or None when matplotlib is not installed"""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        return plt
    except ImportError:
        return None
