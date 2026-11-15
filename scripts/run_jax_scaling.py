import numpy as np
from geocusp.jax_accel import available, two_loop_scalar, loglog_fit
if not available():
    raise SystemExit('JAX not installed; pip install -e .[jax]')
import jax
import jax.numpy as jnp
rho=np.logspace(-1,-5,100)
# Generic positive finite coordinates and generic blow-up directions.
x1,x2,z1,z2=1.2,1.7,.9,1.3
a,b,c,d=.8,1.1,1.5,.6
f=jax.jit(two_loop_scalar)
vals=[]
for r in rho:
    vals.append(float(f(x1,r*a,z1,r*b,x2,r*c,z2,r*d)))
print('JAX devices:',jax.devices())
print('power,prefactor=',loglog_fit(rho,np.asarray(vals)))
