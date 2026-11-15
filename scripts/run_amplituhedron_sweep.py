#!/usr/bin/env python
"""Parameter sweep for the user's amplituhedron pushforward proposal.

Runs several (k,n) cases, two G_+ charts where available, computes positivity
statistics, affine/Pluecker support diagnostics, and saves a machine-readable
report plus figures.  This is a support-geometry experiment, not a claim that
flat Monte Carlo sampling is the exact canonical pushforward measure.
"""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
from geocusp.amplituhedron import (
    user_equation_pushforward, standard_amplituhedron_pushforward,
    pca_coordinates, affine_slice_k1, external_affine_vertices_k1,
)
from geocusp.visualization import (
    positivity_plot, plucker_pca_plot, plucker_pair_plot,
    amplituhedron_affine_plot, boundary_weight_plot,
)


def one_case(k: int, n: int, samples: int, seed: int, out: Path, fig: Path, chart: str):
    r = user_equation_pushforward(k, n, samples, seed, chart=chart, weights="canonical_proxy")
    Q, S, Vt = pca_coordinates(r.plucker)
    outrec = {
        "k": k, "n": n, "samples": samples, "seed": seed, "chart": chart,
        "plucker_dimension": int(r.plucker.shape[1]),
        "min_positivity_margin": float(r.positive_margins.min()),
        "median_positivity_margin": float(np.median(r.positive_margins)),
        "pca_singular_values": [float(x) for x in S[:min(8, len(S))]],
    }
    stem=f"amp_k{k}_n{n}_{chart}"
    positivity_plot(r.positive_margins, fig/f"{stem}_positivity.png",
                    title=fr"User pushforward $G_+({k},{n})$: positivity margin")
    if r.plucker.shape[1] >= 3:
        plucker_pca_plot(r.plucker, fig/f"{stem}_plucker_pca.png",
                         title=fr"User pushforward $G_+({k},{n})$: Plücker PCA")
    if r.plucker.shape[1] >= 2:
        plucker_pair_plot(r.plucker, r.labels, 0, 1, r.positive_margins,
                          fig/f"{stem}_plucker_pair.png",
                          title=fr"$G_+({k},{n})$: Plücker pair / boundary margin")
    if r.weights is not None:
        boundary_weight_plot(r.positive_margins, r.weights,
                              fig/f"{stem}_proxy_weight.png")
    if k == 1 and n == 5:
        pts, hull = affine_slice_k1(r)
        ext = external_affine_vertices_k1(r.Z)
        amplituhedron_affine_plot(pts, hull, ext, fig/f"{stem}_affine_hull.png",
                                  r"User equation pushforward: $k=1,n=5,m=4$")
        outrec["affine_hull_vertices"] = int(len(ext))
    return outrec


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--samples', type=int, default=8000)
    ap.add_argument('--seed', type=int, default=20260918)
    ap.add_argument('--out', default='results')
    ap.add_argument('--figdir', default='figures')
    ap.add_argument('--fast', action='store_true')
    args=ap.parse_args()
    out=Path(args.out); fig=Path(args.figdir); out.mkdir(parents=True,exist_ok=True); fig.mkdir(parents=True,exist_ok=True)
    s=max(1000, args.samples//2) if args.fast else args.samples
    cases=[]
    for idx,(k,n) in enumerate([(1,5),(2,6),(2,7)]):
        cases.append(one_case(k,n,s,args.seed+10*idx,out,fig,'gauge'))
    # Independent chart stress-test for G_+(2,6)
    cases.append(one_case(2,6,s,args.seed+99,out,fig,'polar'))
    rep={"suite":"user_amplituhedron_pushforward_sweep","cases":cases,
         "notes":"Support sampling + positivity diagnostics; canonical_proxy is visualization-only."}
    (out/'amplituhedron_sweep_report.json').write_text(json.dumps(rep,indent=2))
    print(json.dumps(rep,indent=2))

if __name__=='__main__': main()
