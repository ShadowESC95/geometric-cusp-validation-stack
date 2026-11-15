from __future__ import annotations
import numpy as np

def _is_int_line(s):
    try:
        [int(x) for x in s.split()]
        return True
    except Exception:
        return False

def read_locpot(path):
    """Read a standard single-field VASP LOCPOT/CHGCAR volumetric dataset.

    Returns lattice vectors, grid dimensions and the scalar field. PAW/augmentation
    tails after the first grid block are ignored.
    """
    lines=open(path,'r',encoding='utf-8',errors='ignore').read().splitlines()
    i=0
    comment=lines[i]; i+=1
    scale=float(lines[i]); i+=1
    lattice=np.array([[float(x) for x in lines[i+j].split()] for j in range(3)])*scale; i+=3
    if _is_int_line(lines[i].strip()):
        counts=list(map(int,lines[i].split())); i+=1
    else:
        species=lines[i].split(); i+=1
        counts=list(map(int,lines[i].split())); i+=1
    n_atoms=sum(counts)
    if i<len(lines) and lines[i].strip().lower().startswith('selective'):
        i+=1
    i+=1  # Direct / Cartesian
    i+=n_atoms
    while i<len(lines) and not lines[i].strip(): i+=1
    if i>=len(lines): raise ValueError('No volumetric grid dimensions found')
    dims=np.array(list(map(int,lines[i].split())),dtype=int); i+=1
    total=int(np.prod(dims)); vals=[]
    while i<len(lines) and len(vals)<total:
        vals.extend(float(x) for x in lines[i].split()); i+=1
    if len(vals)<total: raise ValueError(f'Expected {total} field values, got {len(vals)}')
    field=np.asarray(vals[:total],dtype=float).reshape(tuple(dims))
    return lattice,dims,field
