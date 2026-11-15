from __future__ import annotations

import itertools
from typing import Iterable

import numpy as np


def ordered_minors(C: np.ndarray, k: int) -> dict[tuple[int, ...], float]:
    """All ordered kxk minors of a single matrix."""
    if C.ndim != 2 or C.shape[0] != k:
        raise ValueError(f"expected shape ({k}, n), got {C.shape}")
    return {I: float(np.linalg.det(C[:, I])) for I in itertools.combinations(range(C.shape[1]), k)}


def positivity_margin(C: np.ndarray, k: int) -> float:
    vals = np.fromiter(ordered_minors(C, k).values(), dtype=float)
    return float(vals.min())


def is_totally_positive(C: np.ndarray, k: int, tol: float = 0.0) -> bool:
    return positivity_margin(C, k) > tol


def sample_gplus_k1(n: int, samples: int, rng: np.random.Generator) -> np.ndarray:
    """Exact gauge-fixed chart for G_+(1,n)/GL(1): positive simplex."""
    if n < 2 or samples < 1:
        raise ValueError("n>=2 and samples>=1 required")
    return rng.dirichlet(np.ones(n), size=samples)


def sample_gplus_k2_gauge(n: int, samples: int, rng: np.random.Generator,
                          log_scale: float = 0.7,
                          slope_spread: float = 1.0) -> np.ndarray:
    """Gauge-fixed full-dimensional positive chart for G_+(2,n).

    Fix C_1=e_1, C_2=e_2. For j>=3 use

        C_j = (-r_j y_j, y_j)^T,

    with y_j>0 and r_3<r_4<...<r_n. Then

        Delta_1j=y_j>0,
        Delta_2j=r_j y_j>0,
        Delta_ij=y_i y_j (r_j-r_i)>0.

    This has exactly 2(n-2)=dim G(2,n) parameters, so no residual GL(2)
    redundancy remains in the sampled chart.
    """
    if n < 3 or samples < 1:
        raise ValueError("n>=3 and samples>=1 required")
    y = np.exp(rng.normal(0.0, log_scale, size=(samples, n - 2)))
    # Positive ordered gaps generate strictly increasing r_j.
    gaps = np.exp(rng.normal(0.0, slope_spread, size=(samples, n - 2)))
    r = np.cumsum(gaps, axis=1)
    C = np.zeros((samples, 2, n), dtype=float)
    C[:, 0, 0] = 1.0
    C[:, 1, 1] = 1.0
    C[:, 0, 2:] = -r * y
    C[:, 1, 2:] = y
    return C


def sample_gplus_k2_polar(n: int, samples: int, rng: np.random.Generator,
                          theta_width: float = 0.92 * np.pi) -> np.ndarray:
    """Alternative positive chart using ordered angles; useful as an independent check.

    Columns are r_i(cos theta_i, sin theta_i) with ordered theta_i in an interval of
    width < pi. The representation has GL(2) redundancy and is therefore best used as
    an independent positivity/support sampler rather than as the primary gauge chart.
    """
    out = np.empty((samples, 2, n), dtype=float)
    for s in range(samples):
        gaps = rng.dirichlet(np.ones(n + 1)) * theta_width
        theta = np.cumsum(gaps[:-1])
        theta -= theta.mean()
        radii = np.exp(rng.normal(0.0, 0.7, n))
        out[s, 0, :] = radii * np.cos(theta)
        out[s, 1, :] = radii * np.sin(theta)
    return out


def sample_gplus_k2(n: int, samples: int, rng: np.random.Generator, *,
                    chart: str = "gauge") -> np.ndarray:
    if chart == "gauge":
        return sample_gplus_k2_gauge(n, samples, rng)
    if chart == "polar":
        return sample_gplus_k2_polar(n, samples, rng)
    raise ValueError(f"unknown chart {chart!r}")


def batch_minors(C: np.ndarray, k: int) -> tuple[np.ndarray, list[tuple[int, ...]]]:
    """Vectorized minor table for a batch of kxn matrices."""
    if C.ndim != 3 or C.shape[1] != k:
        raise ValueError("expected batch shape (samples,k,n)")
    labels = list(itertools.combinations(range(C.shape[2]), k))
    P = np.empty((C.shape[0], len(labels)), dtype=float)
    for j, I in enumerate(labels):
        P[:, j] = np.linalg.det(C[:, :, I])
    return P, labels
