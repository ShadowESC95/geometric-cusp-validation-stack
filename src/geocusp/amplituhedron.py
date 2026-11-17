from __future__ import annotations

import itertools
from dataclasses import dataclass
from typing import Optional

import numpy as np
from scipy.spatial import ConvexHull

from .positive_grassmannian import batch_minors, positivity_margin, sample_gplus_k1, sample_gplus_k2


@dataclass
class PushforwardResult:
    C: np.ndarray
    Z: np.ndarray
    Y: np.ndarray
    positive_margins: np.ndarray
    plucker: np.ndarray
    labels: list[tuple[int, ...]]
    weights: Optional[np.ndarray] = None


@dataclass
class AmplituhedronCloud:
    result: PushforwardResult
    affine_points: Optional[np.ndarray] = None
    hull: Optional[ConvexHull] = None
    metadata: Optional[dict] = None


def positive_external_Z(n: int, d: int, t_min: float = 0.2, t_max: float = 2.4) -> np.ndarray:
    """Moment-curve external data Z_i=(1,t_i,...,t_i^{d-1})."""
    if n < d:
        raise ValueError("n must be >= d")
    t = np.linspace(t_min, t_max, n)
    return np.vstack([t ** p for p in range(d)]).T


def canonical_proxy_weights(C: np.ndarray, k: int, floor: float = 1e-14) -> np.ndarray:
    """Positive-boundary emphasis weight 1/prod(minors).

    This is intentionally named *proxy*: it is not asserted to be the exact pushforward
    density of the canonical form. It is useful for visualizing how positive-boundary
    concentration behaves under a chosen chart.
    """
    P, _ = batch_minors(C, k)
    return 1.0 / np.maximum(np.prod(P, axis=1), floor)


def pminor_matrix(Y: np.ndarray, k: int) -> tuple[np.ndarray, list[tuple[int, ...]]]:
    labels = list(itertools.combinations(range(Y.shape[1]), k))
    p = np.array([np.linalg.det(Y[:, I]) for I in labels], dtype=float)
    scale = np.linalg.norm(p)
    return (p / scale if scale > 0 else p), labels


def pushforward(C: np.ndarray, Z: np.ndarray, k: int, weights: Optional[np.ndarray] = None) -> PushforwardResult:
    Ys = np.einsum("ski,id->skd", C, Z)
    P = []
    labels = None
    for Y in Ys:
        pl, labels = pminor_matrix(Y, k)
        P.append(pl)
    P = np.asarray(P)
    margins = np.asarray([positivity_margin(c, k) for c in C])
    return PushforwardResult(C=C, Z=Z, Y=Ys, positive_margins=margins,
                             plucker=P, labels=labels or [], weights=weights)


def user_equation_pushforward(k: int, n: int, samples: int, seed: int = 20260918,
                              chart: str = "gauge", weights: str = "flat") -> PushforwardResult:
    """Numerically realize the user's formal pushforward support.

    User template:

        Omega^(k,n) ~ int dC/Vol[GL(k)] delta(Y-C Z) I(C).

    The code samples a gauge-fixed positive chart and pushes each C through Y=CZ.
    The support cloud is the primary object. Optional `canonical_proxy` weights are
    provided for boundary-emphasized visualization, but are explicitly not identified
    with the exact canonical-form density.
    """
    rng = np.random.default_rng(seed)
    d = 4
    Z = positive_external_Z(n, d)
    if k == 1:
        C = sample_gplus_k1(n, samples, rng)[:, None, :]
    elif k == 2:
        C = sample_gplus_k2(n, samples, rng, chart=chart)
    else:
        raise NotImplementedError("simulation charts currently implemented for k=1,2")
    W = None if weights == "flat" else canonical_proxy_weights(C, k)
    return pushforward(C, Z, k, weights=W)


def standard_amplituhedron_pushforward(k: int, n: int, samples: int, seed: int = 20260918,
                                       chart: str = "gauge") -> PushforwardResult:
    """Standard m=4 amplituhedron-style pushforward with positive external data."""
    rng = np.random.default_rng(seed)
    d = k + 4
    Z = positive_external_Z(n, d)
    if k == 1:
        C = sample_gplus_k1(n, samples, rng)[:, None, :]
    elif k == 2:
        C = sample_gplus_k2(n, samples, rng, chart=chart)
    else:
        raise NotImplementedError("simulation charts currently implemented for k=1,2")
    return pushforward(C, Z, k)


def affine_slice_k1(result: PushforwardResult) -> tuple[np.ndarray, ConvexHull]:
    if result.Y.shape[1:] != (1, result.Z.shape[1]) or result.Z.shape[1] < 4:
        raise ValueError("requires k=1 and embedding dimension >=4")
    Y = result.Y[:, 0, :]
    X = Y[:, 1:4] / Y[:, [0]]
    hull = ConvexHull(X)
    return X, hull


def external_affine_vertices_k1(Z: np.ndarray) -> np.ndarray:
    if Z.shape[1] < 4:
        raise ValueError("need d>=4")
    return Z[:, 1:4] / Z[:, [0]]


def pca_coordinates(P: np.ndarray, n_components: int = 3) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    X = P - P.mean(axis=0)
    _, S, Vt = np.linalg.svd(X, full_matrices=False)
    n = min(n_components, Vt.shape[0])
    Q = X @ Vt[:n].T
    return Q, S, Vt


def barycentric_reconstruction_residual(result: PushforwardResult) -> np.ndarray:
    """For k=1, verify Y lies in convex hull of external points exactly via C weights."""
    if result.C.shape[1] != 1:
        raise ValueError("k=1 only")
    C = result.C[:, 0, :]
    Y = result.Y[:, 0, :]
    return np.max(np.abs(Y - C @ result.Z), axis=1)


# --------------------------------------------------------------------------
# An independent membership test: sign flips
# --------------------------------------------------------------------------
def _bracket(Y: np.ndarray, Z: np.ndarray, indices) -> float:
    return float(np.linalg.det(np.vstack([Y, Z[list(indices)]])))


def sign_flip_membership(Y: np.ndarray, Z: np.ndarray, k: int) -> dict:
    """Test whether Y (k x (k+4)) lies in the m = 4 amplituhedron of Z without
    using a matrix C: the sign-flip characterisation of Arkani-Hamed, Thomas
    and Trnka (arXiv:1704.05069).

      <Y i i+1 j j+1> > 0 for i < j < n,   (-1)^(k-1) <Y i i+1 n 1> > 0,
      and the sequence <Y 1 2 3 i>, i = 4..n, has exactly k sign flips.

    This is independent of the construction Y = C Z, so it is a real test of
    the sampler, and it rejects points outside the amplituhedron.
    """
    n = Z.shape[0]
    boundary = True
    for i in range(n - 1):
        for j in range(i + 2, n - 1):
            boundary &= _bracket(Y, Z, (i, i + 1, j, j + 1)) > 0
        if 1 <= i <= n - 3:
            boundary &= (-1) ** (k - 1) * _bracket(Y, Z, (i, i + 1, n - 1, 0)) > 0
    seq = [_bracket(Y, Z, (0, 1, 2, i)) for i in range(3, n)]
    flips = sum(1 for a, b in zip(seq, seq[1:]) if a * b < 0)
    return {"boundary_positive": bool(boundary), "flips": flips,
            "inside": bool(boundary and flips == k)}
