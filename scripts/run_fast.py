#!/usr/bin/env python
import subprocess, sys
cmds=[
 [sys.executable,'-m','geocusp.amplituhedron_cli','--mode','both','--samples','2000'],
 [sys.executable,'-m','geocusp.one_loop_cli','--run-all','--contour-points','1024'],
 [sys.executable,'-m','geocusp.two_loop_cli','--mode','all'],
 [sys.executable,'-m','geocusp.continuum_cli','--n','48','--nt','50'],
 [sys.executable,'scripts/run_linking.py'],
]
for c in cmds: subprocess.run(c,check=True)
