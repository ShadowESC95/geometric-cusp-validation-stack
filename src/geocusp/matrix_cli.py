from __future__ import annotations
import argparse, json
from pathlib import Path
from .matrix_model import langevin
from .visualization import matrix_plot

def main():
    p=argparse.ArgumentParser(); p.add_argument('--N',type=int,default=8); p.add_argument('--D',type=int,default=4); p.add_argument('--g',type=float,default=.7); p.add_argument('--steps',type=int,default=4000); p.add_argument('--thermalization',type=int,default=1000); p.add_argument('--out',default='results'); a=p.parse_args()
    Path(a.out).mkdir(exist_ok=True); Path('figures').mkdir(exist_ok=True)
    r=langevin(a.N,a.D,a.g,steps=a.steps,thermalization=a.thermalization)
    matrix_plot(r.steps,r.commutator_density,r.eigenvalues,'figures/matrix_commutator_trajectory.png','figures/matrix_eigenvalues.png')
    report={'N':a.N,'D':a.D,'g':a.g,'mean_commutator_density':float(r.commutator_density.mean()),'final_commutator_density':float(r.commutator_density[-1])}
    (Path(a.out)/'matrix_report.json').write_text(json.dumps(report,indent=2)); print(json.dumps(report,indent=2))
if __name__=='__main__': main()
