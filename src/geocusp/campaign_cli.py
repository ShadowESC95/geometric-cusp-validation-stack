from __future__ import annotations
import argparse, subprocess, sys

def main():
    p=argparse.ArgumentParser(); p.add_argument('--fast',action='store_true'); a=p.parse_args()
    cmds=[['amplituhedron_cli','--mode','both','--samples','3000' if a.fast else '8000'],['one_loop_cli','--contour-points','2048' if a.fast else '4096'],['two_loop_cli','--mode','all'],['continuum_cli','--n','64' if a.fast else '96','--nt','80' if a.fast else '120']]
    for mod in cmds:
        subprocess.run([sys.executable,'-m','geocusp.'+mod[0],*mod[1:]],check=True)
if __name__=='__main__': main()
