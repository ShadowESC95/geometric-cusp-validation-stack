from __future__ import annotations

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt


def savefig(fig, path: str | Path, dpi: int = 180) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def amplituhedron_affine_plot(points, hull, external, path, title):
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection="3d")
    ax.scatter(points[:, 0], points[:, 1], points[:, 2], s=5, alpha=0.15, label="pushforward cloud")
    if hull is not None:
        for tri in hull.simplices:
            idx = np.r_[tri, tri[0]]
            ax.plot(points[idx, 0], points[idx, 1], points[idx, 2], lw=0.7, alpha=0.55)
    if external is not None:
        ax.scatter(external[:, 0], external[:, 1], external[:, 2], s=70, marker="^", label="external Z vertices")
        for i, p in enumerate(external):
            ax.text(*p, f" Z{i+1}", fontsize=8)
    ax.set_xlabel(r"$Y_1/Y_0$")
    ax.set_ylabel(r"$Y_2/Y_0$")
    ax.set_zlabel(r"$Y_3/Y_0$")
    ax.set_title(title)
    ax.legend(loc="best")
    savefig(fig, path)


def plucker_pca_plot(P, path, title, highlight=None):
    X = P - P.mean(axis=0)
    _, _, Vt = np.linalg.svd(X, full_matrices=False)
    Q = X @ Vt[:3].T
    fig = plt.figure(figsize=(9, 8))
    ax = fig.add_subplot(111, projection="3d")
    ax.scatter(Q[:, 0], Q[:, 1], Q[:, 2], s=5, alpha=0.18)
    if highlight is not None:
        ax.scatter(Q[highlight, 0], Q[highlight, 1], Q[highlight, 2], s=14, alpha=0.8, label="near-boundary samples")
        ax.legend(loc="best")
    ax.set_xlabel("Plücker PC1")
    ax.set_ylabel("Plücker PC2")
    ax.set_zlabel("Plücker PC3")
    ax.set_title(title)
    savefig(fig, path)


def plucker_pair_plot(P, labels, i, j, margins, path, title):
    fig, ax = plt.subplots(figsize=(8, 7))
    sc = ax.scatter(P[:, i], P[:, j], c=np.log10(np.maximum(margins, 1e-14)), s=5, alpha=0.3)
    ax.set_xlabel(f"Δ{labels[i]}")
    ax.set_ylabel(f"Δ{labels[j]}")
    ax.set_title(title)
    cb = fig.colorbar(sc, ax=ax)
    cb.set_label(r"$\log_{10}$ positivity margin")
    savefig(fig, path)


def positivity_plot(margins, path, title="G_+ positivity margin"):
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(margins, bins=80)
    ax.set_xlabel("minimum ordered k-minor")
    ax.set_ylabel("count")
    ax.set_title(title)
    savefig(fig, path)


def boundary_weight_plot(margins, weights, path):
    order = np.argsort(margins)
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.loglog(np.maximum(margins[order], 1e-15), np.maximum(weights[order], 1e-15), ".", alpha=0.4)
    ax.set_xlabel("positivity margin")
    ax.set_ylabel("canonical-proxy weight")
    ax.set_title("Boundary concentration diagnostic (proxy weight only)")
    savefig(fig, path)


def residue_map_plot(eta, z, path):
    fig, ax = plt.subplots(figsize=(7, 7))
    ax.plot(z.real, z.imag, lw=1.6)
    ax.scatter([0], [0], s=35)
    ax.set_aspect("equal", "box")
    ax.set_xlabel(r"Re$(z)$")
    ax.set_ylabel(r"Im$(z)$")
    ax.set_title(r"Exceptional-divisor map $z=(z_0/x)\eta$")
    savefig(fig, path)


def linking_plot(z, path, winding):
    fig, ax = plt.subplots(figsize=(7, 7))
    ax.plot(z.real, z.imag, lw=1.6)
    ax.scatter([0], [0], s=40, marker="x")
    ax.set_aspect("equal", "box")
    ax.set_xlabel("Re(z)")
    ax.set_ylabel("Im(z)")
    ax.set_title(f"Linking image in z-plane; winding ≈ {winding:.6f}")
    savefig(fig, path)


def radial_scaling_plot(rho, amp, path, fitted_slope=None, reference=None, title="Collective cusp radial scaling"):
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.loglog(rho, np.abs(amp), "o-", label="measured")
    if fitted_slope is not None:
        c = np.abs(amp[-1]) / rho[-1] ** fitted_slope
        ax.loglog(rho, c * rho ** fitted_slope, "--", label=f"fit slope {fitted_slope:.3f}")
    if reference is not None:
        c = np.abs(amp[-1]) / rho[-1] ** reference
        ax.loglog(rho, c * rho ** reference, ":", label=f"reference slope {reference:g}")
    ax.set_xlabel(r"$\rho$")
    ax.set_ylabel("magnitude")
    ax.set_title(title)
    ax.legend()
    savefig(fig, path)


def contour_integrand_plot(phi, vals, path, title):
    fig, axes = plt.subplots(2, 1, figsize=(9, 8), sharex=True)
    axes[0].plot(phi, vals.real)
    axes[0].set_ylabel("Re")
    axes[1].plot(phi, vals.imag)
    axes[1].set_ylabel("Im")
    axes[1].set_xlabel(r"$\varphi$")
    fig.suptitle(title)
    savefig(fig, path)


def matrix_plot(step, comm, eig, path1, path2):
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(step, comm)
    ax.set_xlabel("measurement step")
    ax.set_ylabel(r"$-\sum_{\mu<\nu}\mathrm{Tr}[X_\mu,X_\nu]^2/N$")
    ax.set_title("IKKT-like matrix Langevin trajectory")
    savefig(fig, path1)
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(eig, bins=60)
    ax.set_xlabel("eigenvalue")
    ax.set_ylabel("count")
    ax.set_title("Matrix eigenvalue spectrum")
    savefig(fig, path2)


def entropy_surface(x, t, S, source, path):
    fig, ax = plt.subplots(figsize=(9, 6))
    im = ax.imshow(S, extent=[x[0], x[-1], t[-1], t[0]], aspect="auto")
    fig.colorbar(im, ax=ax, label=r"$S_{\rm OB}$ proxy")
    ax.set_xlabel("x")
    ax.set_ylabel(r"$\tau$")
    ax.set_title(r"Boundary entropy-transport proxy; $\sigma_{IR}$ source")
    savefig(fig, path)
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(x, source)
    ax.set_xlabel("x")
    ax.set_ylabel(r"$\sigma_{IR}(x)$")
    ax.set_title("IR source profile")
    savefig(fig, Path(path).with_name(Path(path).stem + "_source.png"))


def affine_projection_grid(points, external, path, title="Affine amplituhedron projections"):
    """2D projection panel for the k=1 affine image; useful for boundary shape inspection."""
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    pairs = [(0, 1), (0, 2), (1, 2)]
    labels = [(r"$Y_1/Y_0$", r"$Y_2/Y_0$"), (r"$Y_1/Y_0$", r"$Y_3/Y_0$"), (r"$Y_2/Y_0$", r"$Y_3/Y_0$")]
    for ax, (i,j), lab in zip(axes, pairs, labels):
        ax.hexbin(points[:,i], points[:,j], gridsize=55, mincnt=1)
        if external is not None:
            ax.scatter(external[:,i], external[:,j], s=35, marker='^')
        ax.set_xlabel(lab[0]); ax.set_ylabel(lab[1])
    fig.suptitle(title)
    savefig(fig, path)


def alpha_z_path_plot(alpha, z, path):
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(alpha.real, alpha.imag, label=r"$\alpha$")
    ax.plot(z.real, z.imag, label=r"$z(\alpha)$")
    ax.scatter([0], [0], marker='x')
    ax.set_aspect('equal', 'box')
    ax.set_xlabel('real'); ax.set_ylabel('imaginary')
    ax.set_title(r"Strict-cusp modulus versus induced cross-ratio path")
    ax.legend()
    savefig(fig, path)


def bracket_scaling_plot(rho, curves, path):
    fig, ax = plt.subplots(figsize=(8, 5))
    for label, vals in curves.items():
        ax.loglog(rho, np.abs(vals), 'o-', label=label)
    ax.set_xlabel(r"$\rho$"); ax.set_ylabel('absolute bracket / scalar'); ax.set_title('Singular-sector scaling diagnostics'); ax.legend()
    savefig(fig, path)
