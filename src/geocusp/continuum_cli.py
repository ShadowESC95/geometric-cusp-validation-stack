from __future__ import annotations
import argparse, json
from pathlib import Path
from .continuum import run_entropy_transport
from .visualization import entropy_surface

def main():
    p=argparse.ArgumentParser(); p.add_argument('--n',type=int,default=96); p.add_argument('--nt',type=int,default=120); p.add_argument('--out',default='results'); a=p.parse_args(); Path(a.out).mkdir(exist_ok=True); Path('figures').mkdir(exist_ok=True)
    r=run_entropy_transport(n=a.n,nt=a.nt)
    entropy_surface(r.x,r.t,r.S,r.sigma,'figures/entropy_transport_proxy.png')
    report={'max_S':float(r.S.max()),'final_integrated_S':float(r.S[-1].sum()*(r.x[1]-r.x[0])),'max_sigma':float(r.sigma.max())}; (Path(a.out)/'continuum_report.json').write_text(json.dumps(report,indent=2)); print(json.dumps(report,indent=2))
if __name__=='__main__': main()
