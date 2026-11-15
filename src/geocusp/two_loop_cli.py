from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from .negative_geometry import symbolic_two_loop_integrand, two_loop_blowup_leading_coefficient, collective_scaling, iterated_residue_order_scan
from .visualization import radial_scaling_plot


def numerical_two_loop_scaling(rho_values, sample=1.0):
    """Evaluate the exact scalar integrand on a fixed positive blow-up ray.

    This is an intentionally simple diagnostic, not a regulated two-loop integral.
    """
    import sympy as sp
    vars_, I, _, _ = symbolic_two_loop_integrand()
    x1, y1, z1, w1, x2, y2, z2, w2 = vars_
    subs_fixed = {x1: 1.1, z1: 0.8, x2: 1.7, z2: 1.3,
                  y1: sp.Symbol('rho') * 0.6, w1: sp.Symbol('rho') * 1.2,
                  y2: sp.Symbol('rho') * 0.9, w2: sp.Symbol('rho') * 1.5}
    f = sp.lambdify(sp.Symbol('rho'), I.subs(subs_fixed), "numpy")
    return np.asarray(rho_values, dtype=float), np.asarray([float(f(r)) for r in rho_values])


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--mode", choices=["exact", "scaling", "all"], default="all")
    p.add_argument("--out", default="results")
    p.add_argument("--figdir", default="figures")
    a = p.parse_args()
    Path(a.out).mkdir(parents=True, exist_ok=True)
    Path(a.figdir).mkdir(parents=True, exist_ok=True)
    report = {}
    if a.mode in ("exact", "all"):
        vars_, I, X, n = symbolic_two_loop_integrand()
        scaled, lead = two_loop_blowup_leading_coefficient()
        report["two_loop_integrand"] = str(I)
        report["mutual_bracket"] = str(X)
        report["numerator"] = str(n)
        report["rho4_scalar_lead"] = str(lead)
    if a.mode in ("scaling", "all"):
        data = []
        for L in range(1, 8):
            graph = [(i, i + 1) for i in range(L - 1)]
            data.append(collective_scaling(L, graph))
        report["tree_scaling"] = data
        report["residue_order_scan"] = iterated_residue_order_scan()
        rho = np.logspace(-2, -6, 17)
        rr, vals = numerical_two_loop_scaling(rho)
        slope = float(np.polyfit(np.log(rr), np.log(np.abs(vals)), 1)[0])
        report["fixed_ray_numeric_slope"] = slope
        radial_scaling_plot(rr, vals, Path(a.figdir) / "07_two_loop_fixed_ray_scaling.png",
                            fitted_slope=slope, reference=-4.0,
                            title=r"Exact two-loop scalar blow-up: $I_4^{(2)}(\rho)$")
    (Path(a.out) / "two_loop_audit.json").write_text(json.dumps(report, indent=2, default=str))
    print(json.dumps(report, indent=2, default=str))


if __name__ == "__main__":
    main()
