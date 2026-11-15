from __future__ import annotations
import numpy as np

def available():
    try:
        import jax, jax.numpy as jnp
        return True
    except Exception:
        return False


def two_loop_scalar(xs, ys, zs, ws, x2s, y2s, z2s, w2s):
    import jax.numpy as jnp
    n=ws*xs+w2s*x2s+ys*zs+y2s*z2s
    D1=ws*xs*ys*zs; D2=w2s*x2s*y2s*z2s
    X=-(xs-x2s)*(ws-w2s)-(ys-y2s)*(zs-z2s)
    return -n/(D1*D2*X)


def loglog_fit(rho, y):
    p=np.polyfit(np.log(rho),np.log(np.abs(y)),1)
    return float(p[0]),float(np.exp(p[1]))
