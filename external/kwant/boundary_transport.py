"""Kwant transport probe for a spatially varying boundary mass."""
from __future__ import annotations
import numpy as np


def make_system(L=80, W=12, hopping=1.0, mass0=0.4, edge_amp=1.2):
    try:
        import kwant
    except Exception as exc:
        raise RuntimeError("Kwant unavailable") from exc
    lat=kwant.lattice.square(norbs=2)
    syst=kwant.Builder()
    sigma_z=np.array([[1,0],[0,-1]],complex)
    sigma_x=np.array([[0,1],[1,0]],complex)
    def onsite(site):
        x,y=site.pos
        m=mass0+edge_amp*np.exp(-(x-L/2)**2/(0.12*L)**2)
        return m*sigma_z
    syst[lat.shape((lambda p: 0<=p[0]<L and 0<=p[1]<W),(0,0))]=onsite
    syst[lat.neighbors()] = -hopping*np.eye(2)
    lead=kwant.Builder(kwant.TranslationalSymmetry((-1,0)))
    lead[lat.shape(lambda p: 0<=p[1]<W,(0,0))]=lambda site: mass0*sigma_z
    lead[lat.neighbors()] = -hopping*np.eye(2)
    syst.attach_lead(lead); syst.attach_lead(lead.reversed())
    return syst.finalized()


def transmission(E=0.2):
    import kwant
    s=make_system()
    smat=kwant.smatrix(s,E)
    return float(smat.transmission(1,0))

if __name__=='__main__': print('T(E=0.2)=', transmission())
