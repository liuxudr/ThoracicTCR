"""Repertoire diversity & clonality metrics."""

from __future__ import annotations

import numpy as np
import pandas as pd


def _counts(clones: pd.DataFrame, count_col: str = "duplicate_count") -> np.ndarray:
    c = clones[count_col].to_numpy(dtype=float)
    return c[c > 0]


def _proportions(counts: np.ndarray) -> np.ndarray:
    total = counts.sum()
    return counts / total if total > 0 else counts


def shannon_entropy(counts: np.ndarray) -> float:
    p = _proportions(counts)
    p = p[p > 0]
    return float(-np.sum(p * np.log(p))) if p.size else 0.0


def simpson_index(counts: np.ndarray) -> float:
    """1 - Σp²  (probability that two reads are different clones)."""
    p = _proportions(counts)
    return float(1.0 - np.sum(p ** 2))


def inverse_simpson(counts: np.ndarray) -> float:
    p = _proportions(counts)
    denom = float(np.sum(p ** 2))
    return float(1.0 / denom) if denom > 0 else 0.0


def pielou_evenness(counts: np.ndarray) -> float:
    """Shannon / log(richness)."""
    h = shannon_entropy(counts)
    s = (counts > 0).sum()
    return float(h / np.log(s)) if s > 1 else 0.0


def gini_coefficient(counts: np.ndarray) -> float:
    if counts.size == 0:
        return 0.0
    x = np.sort(counts)
    n = x.size
    cum = np.cumsum(x)
    return float((n + 1 - 2 * cum.sum() / cum[-1]) / n) if cum[-1] > 0 else 0.0


def chao1(counts: np.ndarray) -> float:
    """Chao1 richness estimator. Requires integer counts."""
    c = counts.astype(int)
    f1 = (c == 1).sum()
    f2 = (c == 2).sum()
    s_obs = (c > 0).sum()
    if f2 > 0:
        return float(s_obs + (f1 * (f1 - 1)) / (2 * (f2 + 1)))
    return float(s_obs + f1 * (f1 - 1) / 2)


def hill_number(counts: np.ndarray, q: float) -> float:
    """Generalized Hill diversity of order q (effective # species)."""
    p = _proportions(counts)
    p = p[p > 0]
    if p.size == 0:
        return 0.0
    if q == 1.0:
        return float(np.exp(-np.sum(p * np.log(p))))
    return float(np.sum(p ** q) ** (1.0 / (1.0 - q)))


def top_n_proportion(counts: np.ndarray, n: int) -> float:
    if counts.size == 0:
        return 0.0
    s = np.sort(counts)[::-1]
    return float(s[:n].sum() / counts.sum())


def d50(counts: np.ndarray) -> int:
    """Number of unique clones required to cover ≥ 50% of repertoire."""
    if counts.size == 0:
        return 0
    s = np.sort(counts)[::-1]
    cum = np.cumsum(s) / s.sum()
    return int(np.searchsorted(cum, 0.5) + 1)


def clonality(counts: np.ndarray) -> float:
    """1 - normalized Shannon ∈ [0,1]; higher = more clonal."""
    h = shannon_entropy(counts)
    s = (counts > 0).sum()
    return float(1.0 - h / np.log(s)) if s > 1 else 0.0


def compute_diversity(
    clones: pd.DataFrame,
    count_col: str = "duplicate_count",
) -> dict:
    """Compute all standard diversity / clonality metrics."""
    c = _counts(clones, count_col)
    n_clones = int((c > 0).sum())
    total_reads = int(c.sum())
    return {
        "n_clones": n_clones,
        "total_reads": total_reads,
        "shannon": shannon_entropy(c),
        "simpson": simpson_index(c),
        "inverse_simpson": inverse_simpson(c),
        "pielou_evenness": pielou_evenness(c),
        "gini": gini_coefficient(c),
        "chao1": chao1(c),
        "hill_q0": float(n_clones),
        "hill_q1": hill_number(c, 1.0),
        "hill_q2": hill_number(c, 2.0),
        "top1_frac": top_n_proportion(c, 1),
        "top10_frac": top_n_proportion(c, 10),
        "top100_frac": top_n_proportion(c, 100),
        "d50": d50(c),
        "clonality": clonality(c),
    }
