from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from .cusp_audit import run_audit
from .residue import blowup_intersection_yw, strict_cusp_blowup_zero_residue, chart_transition_source_to_U12
from .visualization import residue_map_plot, linking_plot, contour_integrand_plot, alpha_z_path_plot
from .contours import contour_points, winding_number, benchmark_F1


def main() -> None:
    p=argparse.ArgumentParser(description="Exact one-loop Grassmannian/cusp audit")
    p.add_argument('--contour-points',type=int,default=4096)
    p.add_argument('--out',default='results')
    p.add_argument('--figdir',default='figures')
    a=p.parse_args()
    Path(a.out).mkdir(parents=True,exist_ok=True); Path(a.figdir).mkdir(parents=True,exist_ok=True)

    res, chi = blowup_intersection_yw()
    pull, zero = strict_cusp_blowup_zero_residue()
    mapping, inverse, J = chart_transition_source_to_U12()

    eta = np.exp(2j*np.pi*np.arange(a.contour_points)/a.contour_points)
    z_ex = 0.9*eta
    residue_map_plot(eta,z_ex,Path(a.figdir)/'one_loop_exceptional_to_z.png')
    linking_plot(z_ex,Path(a.figdir)/'one_loop_linking_z_plane.png',winding_number(z_ex))

    z,phi=contour_points(1.0,a.contour_points)
    vals=np.asarray([benchmark_F1(zz) for zz in z])
    contour_integrand_plot(phi,vals,Path(a.figdir)/'one_loop_F1_contour.png',
                           title=r'$F^{(1)}(z)=\log^2 z+\pi^2$ on $|z|=1$')

    audit=run_audit(n=a.contour_points)
    alpha_z_path_plot(audit.alpha_curve, audit.z_curve, Path(a.figdir)/'one_loop_alpha_vs_z.png')
    report={
      'source_chart_form':'dlog(x) wedge dlog(y) wedge dlog(z0) wedge dlog(w)',
      'exceptional_residue':str(res),
      'exceptional_cross_ratio':str(chi),
      'strict_cusp_blowup_pullback_coefficient':str(pull),
      'strict_cusp_radial_residue':str(zero),
      'source_to_U12_mapping':{str(k):str(v) for k,v in mapping.items()},
      'U12_to_source_inverse':{str(k):str(v) for k,v in inverse.items()},
      'source_to_U12_log_jacobian':str(J),
      'z_linking_winding':float(audit.linking_winding),
      'z_linking_winding_nearest_integer':int(round(audit.linking_winding)),
      'status':{
        'simple_radial_residue_on_strict_cusp_chart':'ZERO',
        'exceptional_divisor_residue_on_adjacent_propagator_intersection':'NONZERO',
        'exceptional_eta_to_cross_ratio':'z=(z0/x)*eta',
        'literal_alpha_equals_z':'NOT_ESTABLISHED',
      }
    }
    (Path(a.out)/'one_loop_audit.json').write_text(json.dumps(report,indent=2,default=str))
    print(json.dumps(report,indent=2,default=str))

if __name__=='__main__': main()
