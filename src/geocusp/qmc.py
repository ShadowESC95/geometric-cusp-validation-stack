from __future__ import annotations
import numpy as np
from scipy.stats import qmc


def randomized_sobol_integral(f, dim, bounds, samples=2**16, replicates=16, seed=20260918):
    bounds=np.asarray(bounds,float)
    vals=[]
    volume=float(np.prod(bounds[:,1]-bounds[:,0]))
    for r in range(replicates):
        eng=qmc.Sobol(d=dim,scramble=True,seed=seed+r)
        u=eng.random_base2(int(np.log2(samples)))
        x=bounds[:,0]+u*(bounds[:,1]-bounds[:,0])
        fx=np.asarray([f(xx) for xx in x],float)
        vals.append(volume*fx.mean())
    vals=np.asarray(vals)
    return float(vals.mean()), float(vals.std(ddof=1)/np.sqrt(len(vals))), vals
