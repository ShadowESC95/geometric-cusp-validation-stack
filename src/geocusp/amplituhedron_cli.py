from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from .amplituhedron import (
    affine_slice_k1,
    barycentric_reconstruction_residual,
    external_affine_vertices_k1,
    pca_coordinates,
    standard_amplituhedron_pushforward,
    user_equation_pushforward,
)
from .visualization import (
    amplituhedron_affine_plot,
    boundary_weight_plot,
    plucker_pair_plot,
    plucker_pca_plot,
    positivity_plot,
)


def run(mode: str, samples: int, seed: int, out: str, figdir: str, chart: str) -> dict:
    outp = Path(out)
    fig = Path(figdir)
    outp.mkdir(parents=True, exist_ok=True)
    fig.mkdir(parents=True, exist_ok=True)
    report: dict = {}

    if mode in ("user", "both"):
        r = user_equation_pushforward(2, 6, samples, seed, chart=chart, weights="canonical_proxy")
        Q, S, Vt = pca_coordinates(r.plucker)
        near = r.positive_margins <= np.quantile(r.positive_margins, 0.05)
        report["user_k2_n6"] = {
            "samples": samples,
            "chart": chart,
            "min_margin": float(r.positive_margins.min()),
            "mean_margin": float(r.positive_margins.mean()),
            "plucker_dim": int(r.plucker.shape[1]),
            "near_boundary_count": int(near.sum()),
            "pca_explained_singular_values": [float(x) for x in S[:3]],
        }
        positivity_plot(r.positive_margins, fig / "01_user_Gplus_positivity.png",
                        title=r"User pushforward: $G_+(2,6)$ positivity margin")
        plucker_pca_plot(r.plucker, fig / "02_user_plucker_pca.png",
                         title=r"User pushforward: Pl" + "ücker PCA", highlight=np.where(near)[0])
        # Choose the first two minors for a boundary-vs-minor diagnostic.
        plucker_pair_plot(r.plucker, r.labels, 0, 1, r.positive_margins,
                          fig / "03_user_plucker_pair.png",
                          title="Two Plücker coordinates coloured by positivity margin")
        boundary_weight_plot(r.positive_margins, r.weights, fig / "04_user_boundary_weight.png")

        r1 = user_equation_pushforward(1, 5, max(5000, samples // 2), seed + 1)
        pts, hull = affine_slice_k1(r1)
        ext = external_affine_vertices_k1(r1.Z)
        residual = barycentric_reconstruction_residual(r1)
        report["user_k1_n5"] = {
            "samples": int(r1.C.shape[0]),
            "max_affine_reconstruction_residual": float(residual.max()),
            "external_vertices": ext.tolist(),
            "hull_simplices": int(len(hull.simplices)),
        }
        amplituhedron_affine_plot(pts, hull, ext, fig / "05_user_A1n4_affine_hull.png",
                                  title=r"User equation pushforward: $k=1,n=5,m=4$ affine slice")

    if mode in ("standard", "both"):
        r = standard_amplituhedron_pushforward(2, 6, samples, seed + 11, chart=chart)
        Q, S, _ = pca_coordinates(r.plucker)
        report["standard_A_6_2_4"] = {
            "samples": samples,
            "chart": chart,
            "min_margin": float(r.positive_margins.min()),
            "plucker_dim": int(r.plucker.shape[1]),
            "singular_values": [float(x) for x in S[:5]],
        }
        plucker_pca_plot(r.plucker, fig / "06_standard_A624_plucker_pca.png",
                         title=r"Standard-style $A_{6,2,4}$ pushforward: Pl" + "ücker PCA")

    (outp / "amplituhedron_report.json").write_text(json.dumps(report, indent=2))
    return report


def main() -> None:
    p = argparse.ArgumentParser(description="Simulate and visualize user-defined and standard-style amplituhedron pushforwards.")
    p.add_argument("--mode", choices=["user", "standard", "both"], default="both")
    p.add_argument("--samples", type=int, default=12000)
    p.add_argument("--seed", type=int, default=20260918)
    p.add_argument("--chart", choices=["gauge", "polar"], default="gauge")
    p.add_argument("--out", default="results")
    p.add_argument("--figdir", default="figures")
    a = p.parse_args()
    report = run(a.mode, a.samples, a.seed, a.out, a.figdir, a.chart)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
