from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import json

import numpy as np

from .amplituhedron import (
    user_equation_pushforward,
    standard_amplituhedron_pushforward,
    affine_slice_k1,
    external_affine_vertices_k1,
    pca_coordinates,
    canonical_proxy_weights,
)
from .visualization import (
    amplituhedron_affine_plot,
    affine_projection_grid,
    plucker_pca_plot,
    positivity_plot,
    plucker_pair_plot,
    boundary_weight_plot,
)


@dataclass
class AmplituhedronStudy:
    user_k1: object
    user_k2: object
    standard_k2: object
    report: dict


def run_study(samples_k2: int = 20000, samples_k1: int = 12000,
              seed: int = 20260918, outdir: str = "results",
              figdir: str = "figures") -> AmplituhedronStudy:
    out = Path(outdir); fig = Path(figdir)
    out.mkdir(parents=True, exist_ok=True); fig.mkdir(parents=True, exist_ok=True)

    r1 = user_equation_pushforward(1, 5, samples_k1, seed)
    pts, hull = affine_slice_k1(r1)
    ext = external_affine_vertices_k1(r1.Z)
    amplituhedron_affine_plot(
        pts, hull, ext, fig / "amplituhedron_user_k1_hull.png",
        r"User pushforward $\Omega^{(1,5)}$: affine amplituhedron hull"
    )
    affine_projection_grid(pts, ext, fig / "amplituhedron_user_k1_projections.png")

    r2 = user_equation_pushforward(2, 6, samples_k2, seed + 1, chart="gauge", weights="canonical_proxy")
    Q2, S2, _ = pca_coordinates(r2.plucker)
    near = r2.positive_margins <= np.quantile(r2.positive_margins, 0.03)
    plucker_pca_plot(r2.plucker, fig / "amplituhedron_user_k2_plucker_pca.png",
                     r"User pushforward $\Omega^{(2,6)}$: Pl" + "ücker PCA", highlight=np.where(near)[0])
    positivity_plot(r2.positive_margins, fig / "amplituhedron_user_k2_positivity.png",
                    title=r"$G_+(2,6)$ positivity margin")
    plucker_pair_plot(r2.plucker, r2.labels, 0, 1, r2.positive_margins,
                      fig / "amplituhedron_user_k2_pair.png",
                      "Two Plücker coordinates vs. boundary proximity")
    if r2.weights is not None:
        boundary_weight_plot(r2.positive_margins, r2.weights,
                             fig / "amplituhedron_user_k2_proxy_weight.png")

    rs = standard_amplituhedron_pushforward(2, 6, samples_k2, seed + 2, chart="gauge")
    _, Ss, _ = pca_coordinates(rs.plucker)
    plucker_pca_plot(rs.plucker, fig / "amplituhedron_standard_A624_pca.png",
                     r"Standard-style $A_{6,2,4}$ pushforward: Pl" + "ücker PCA")

    residual = np.max(np.abs(r1.Y[:, 0, :] - np.einsum("si,id->sd", r1.C[:, 0, :], r1.Z)), axis=1)
    report = {
        "user_formula": r"Omega^(k,n) = int dC/Vol[GL(k)] delta(Y-CZ) I(C)",
        "user_k1_n5": {
            "samples": samples_k1,
            "max_pushforward_residual": float(residual.max()),
            "n_affine_hull_vertices": int(len(ext)),
        },
        "user_k2_n6": {
            "samples": samples_k2,
            "min_positive_margin": float(r2.positive_margins.min()),
            "median_positive_margin": float(np.median(r2.positive_margins)),
            "n_plucker_coordinates": int(r2.plucker.shape[1]),
            "near_boundary_count": int(near.sum()),
            "pca_singular_values": [float(v) for v in S2[:5]],
        },
        "standard_A624": {
            "samples": samples_k2,
            "n_plucker_coordinates": int(rs.plucker.shape[1]),
            "min_positive_margin": float(rs.positive_margins.min()),
            "pca_singular_values": [float(v) for v in Ss[:5]],
        },
    }
    (out / "amplituhedron_deep_report.json").write_text(json.dumps(report, indent=2))
    return AmplituhedronStudy(r1, r2, rs, report)
