from __future__ import annotations

from dataclasses import dataclass
import numpy as np
import sympy as sp

from .residue import blowup_intersection_yw, strict_cusp_blowup_zero_residue, chart_transition_source_to_U12
from .twistor import standard_external_Z4
from .contours import winding_number


@dataclass
class CuspAudit:
    exceptional_residue: str
    exceptional_cross_ratio: str
    strict_radial_residue: str
    transition_jacobian: str
    linking_winding: float
    alpha_curve: np.ndarray
    z_curve: np.ndarray


def alpha_path_z(alpha_values: np.ndarray, rho: float = 1e-5) -> np.ndarray:
    """Evaluate the four-point cross ratio along a smooth off-boundary path approaching the strict line."""
    Z = standard_external_Z4().astype(complex)
    out = []
    for alpha in alpha_values:
        A = Z[1] + rho * Z[3]
        B = Z[0] + alpha * Z[2] + rho * Z[3]
        vals = []
        for i, j in [(0,1),(1,2),(2,3),(0,3)]:
            vals.append(np.linalg.det(np.column_stack([A, B, Z[i], Z[j]])))
        b12,b23,b34,b14 = vals
        out.append(b12*b34/(b23*b14))
    return np.asarray(out)


def run_audit(n: int = 2048, radius: float = 0.9) -> CuspAudit:
    r, chi = blowup_intersection_yw()
    _, zero = strict_cusp_blowup_zero_residue()
    _, _, J = chart_transition_source_to_U12()
    eta = radius * np.exp(2j * np.pi * np.arange(n) / n)
    # Generic nonzero scale factor does not change winding.
    z = (0.9 / 1.3) * eta
    alpha = 0.5 + 0.3 * np.exp(2j * np.pi * np.arange(n) / n)
    z_alpha = alpha_path_z(alpha)
    return CuspAudit(str(r), str(chi), str(zero), str(J), winding_number(z), alpha, z_alpha)
