from __future__ import annotations
import csv
from pathlib import Path
import numpy as np


def export_1d_field(x, fields: dict[str,np.ndarray], path='results/geometric_field.csv'):
    p=Path(path); p.parent.mkdir(parents=True,exist_ok=True)
    names=['x']+list(fields)
    with p.open('w',newline='') as f:
        w=csv.writer(f); w.writerow(names)
        for i,xx in enumerate(x): w.writerow([float(xx)]+[float(fields[n][i]) for n in fields])
    return p


def export_openfoam_table(x,sigma,path='external/openfoam/generated/sigmaIR.xy'):
    p=Path(path); p.parent.mkdir(parents=True,exist_ok=True); p.write_text('\n'.join(f'{xx:.16e} {ss:.16e}' for xx,ss in zip(x,sigma))); return p
