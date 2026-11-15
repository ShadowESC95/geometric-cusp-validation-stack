#!/usr/bin/env python
from __future__ import annotations
import json, subprocess, sys
from pathlib import Path
from geocusp.backends import detect

ROOT=Path(__file__).resolve().parents[1]
D=detect(); report={}

# Python-native backends are run when installed. Heavy/system backends receive instructions.
if D['qutip']['available']:
    r=subprocess.run([sys.executable,'external/qutip/open_system_entropy.py'],cwd=ROOT,text=True,capture_output=True)
    report['qutip']={'returncode':r.returncode,'stdout':r.stdout,'stderr':r.stderr}
if D['fipy']['available']:
    r=subprocess.run([sys.executable,'external/fipy/boundary_entropy.py'],cwd=ROOT,text=True,capture_output=True)
    report['fipy']={'returncode':r.returncode,'stdout':r.stdout,'stderr':r.stderr}
if D['kwant']['available']:
    r=subprocess.run([sys.executable,'-c','from external.kwant.boundary_transport import transmission; print(transmission())'],cwd=ROOT,text=True,capture_output=True)
    report['kwant']={'returncode':r.returncode,'stdout':r.stdout,'stderr':r.stderr}
if D['meep']['available']:
    r=subprocess.run([sys.executable,'external/meep/geometric_boundary_wave.py'],cwd=ROOT,text=True,capture_output=True)
    report['meep']={'returncode':r.returncode,'stdout':r.stdout[-2000:],'stderr':r.stderr[-2000:]}
for name in ['fenicsx','firedrake','deal_ii','comsol','form','vasp','openfoam','sage','macaulay2','singular','cadabra']:
    report[name]={'available':D[name]['available'],'note':'Adapter present; execute native backend with external/ script or case files.'}
Path('results').mkdir(exist_ok=True)
Path('results/optional_backend_runs.json').write_text(json.dumps(report,indent=2))
print(json.dumps(report,indent=2))
