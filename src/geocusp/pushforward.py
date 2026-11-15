from __future__ import annotations
import numpy as np


def gaussian_delta_weight(Y: np.ndarray, CZ: np.ndarray, sigma: float) -> float:
    r=Y-CZ
    k=r.size
    return float(np.exp(-0.5*np.sum(np.abs(r)**2)/sigma**2)/(2*np.pi*sigma**2)**(k/2))


def kde_eval(samples_CZ: np.ndarray, Y: np.ndarray, sigma: float) -> float:
    r=samples_CZ-Y.reshape(1,*Y.shape)
    k=Y.size
    weights=np.exp(-0.5*np.sum(np.abs(r)**2,axis=tuple(range(1,r.ndim)))/sigma**2)
    return float(weights.mean()/(2*np.pi*sigma**2)**(k/2))
