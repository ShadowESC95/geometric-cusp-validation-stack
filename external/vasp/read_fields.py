"""VASP volumetric-field importer for LOCPOT/CHGCAR-derived vacuum/VEV studies."""
from __future__ import annotations
from pathlib import Path
import numpy as np


def read_volumetric(path: str):
    lines=Path(path).read_text(encoding='utf-8',errors='ignore').splitlines()
    i=0; comment=lines[i]; i+=1; scale=float(lines[i]); i+=1
    lattice=np.array([[float(x) for x in lines[i+j].split()] for j in range(3)])*scale; i+=3
    # VASP may list species then counts, or counts directly.
    toks=lines[i].split()
    if all(tok.lstrip('+-').isdigit() for tok in toks): counts=list(map(int,toks)); i+=1
    else:
        i+=1; counts=list(map(int,lines[i].split())); i+=1
    n_atoms=sum(counts)
    if i<len(lines) and lines[i].lower().startswith('selective'): i+=1
    i+=1; i+=n_atoms
    while i<len(lines) and not lines[i].strip(): i+=1
    dims=tuple(map(int,lines[i].split())); i+=1
    total=int(np.prod(dims)); vals=[]
    while i<len(lines) and len(vals)<total: vals.extend(float(x) for x in lines[i].split()); i+=1
    if len(vals)<total: raise ValueError('incomplete volumetric field')
    return lattice,dims,np.asarray(vals[:total]).reshape(dims)
