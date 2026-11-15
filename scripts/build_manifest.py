#!/usr/bin/env python
"""Build a reproducibility manifest with hashes and runtime metadata."""
from __future__ import annotations
from pathlib import Path
import hashlib, json, platform, sys, time, subprocess

ROOT=Path(__file__).resolve().parents[1]
IGNORE={'.pytest_cache','__pycache__','.git'}

def sha256(p: Path):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024), b''): h.update(chunk)
    return h.hexdigest()

def main():
    files=[]
    for p in ROOT.rglob('*'):
        if not p.is_file() or any(part in IGNORE for part in p.parts): continue
        if p.name=='manifest.json': continue
        files.append({'path':str(p.relative_to(ROOT)),'sha256':sha256(p),'bytes':p.stat().st_size})
    try:
        git=subprocess.run(['git','rev-parse','HEAD'],cwd=ROOT,text=True,capture_output=True).stdout.strip()
    except Exception: git=''
    env={k:__import__('importlib').util.find_spec(k) is not None for k in
         ['sympy','numpy','scipy','mpmath','matplotlib','jax','plotly','flint','sageall','petsc4py','dolfinx','firedrake','fipy','meep','kwant','qutip']}
    manifest={'generated_utc':time.time(),'python':sys.version,'platform':platform.platform(),'git_head':git,'modules_available':env,'files':sorted(files,key=lambda x:x['path'])}
    out=ROOT/'results'/'manifest.json'; out.parent.mkdir(exist_ok=True); out.write_text(json.dumps(manifest,indent=2)); print(out)
if __name__=='__main__': main()
