"""Meep FDTD probe driven by an effective boundary profile.

The profile epsilon(x) is a placeholder map from the user's emergent field layer
into electrodynamics; it must be replaced by the derived constitutive relation
before physical claims are made.
"""
from __future__ import annotations
import math


def run(resolution=80, sx=20.0, frequency=0.20, chi0=2.0, width=2.0):
    try:
        import meep as mp
    except Exception as exc:
        raise RuntimeError("Meep unavailable") from exc
    cell=mp.Vector3(sx,0,0)
    pml_layers=[mp.PML(1.0)]
    def eps_func(p):
        return 1.0 + chi0*math.exp(-(p.x/width)**2)
    medium=mp.Medium(epsilon_func=eps_func)
    sources=[mp.Source(mp.GaussianSource(frequency,fwidth=0.05), component=mp.Ey, center=mp.Vector3(-6,0,0), size=mp.Vector3(0,0,0))]
    sim=mp.Simulation(cell_size=cell,boundary_layers=pml_layers,geometry=[mp.Block(mp.Vector3(mp.inf,mp.inf,mp.inf),material=medium)],sources=sources,resolution=resolution)
    sim.run(until=100)
    return sim

if __name__=='__main__': run(); print('Meep run complete')
