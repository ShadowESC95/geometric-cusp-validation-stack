#!/usr/bin/env python
"""One-command reproducible research campaign.

This is the main orchestration layer for the multi-backend validation stack.
Optional scientific backends are detected and reported; they are never silently
replaced by toy models.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]


def run(cmd):
    print('\n>>>', ' '.join(map(str,cmd)), flush=True)
    subprocess.run(cmd, cwd=ROOT, check=True)


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--fast',action='store_true')
    p.add_argument('--with-optional',action='store_true', help='run adapters that are actually installed')
    p.add_argument('--samples-k2',type=int,default=None)
    args=p.parse_args()
    k2=args.samples_k2 or (6000 if args.fast else 20000)
    k1=4000 if args.fast else 12000
    contour=1024 if args.fast else 4096

    run([sys.executable,'scripts/backend_report.py'])
    run([sys.executable,'scripts/run_amplituhedron_deep.py','--samples-k2',str(k2),'--samples-k1',str(k1)])
    run([sys.executable,'-m','geocusp.one_loop_cli','--contour-points',str(contour)])
    run([sys.executable,'-m','geocusp.two_loop_cli','--mode','all'])
    run([sys.executable,'scripts/run_benchmark.py'])
    run([sys.executable,'scripts/run_linking.py'])
    run([sys.executable,'scripts/run_cusp_paths.py'])
    run([sys.executable,'scripts/run_jax_scaling.py'])
    run([sys.executable,'-m','geocusp.matrix_cli','--N','8','--D','4','--g','0.7','--steps',('1200' if args.fast else '5000'),'--thermalization','500'])
    run([sys.executable,'-m','geocusp.continuum_cli','--n',('64' if args.fast else '128'),'--nt',('80' if args.fast else '160')])
    run([sys.executable,'scripts/export_continuum_fields.py'])
    if args.with_optional:
        run([sys.executable,'scripts/run_optional_backends.py'])

    print('\nCampaign complete. Read results/ and figures/.')

if __name__=='__main__': main()
