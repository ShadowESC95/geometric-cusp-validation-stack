"""Compare the two lateral CD parameter chains at the capped-cycle seam."""

from cmath import exp, log as complex_log
from math import log, pi

import numpy as np
import sympy as sp
from scipy.integrate import quad

from geocusp.contours import winding_number

from .complex_cross_ratio import check_projective_pullback, integrated_projective


def certify_moving_pole():
    check_projective_pullback()
    x, t, z, s, v = sp.symbols("x t z s v")
    # y=x*t. The x in the Jacobian cancels the x in z*x+y.
    reduced = 2*(1+z)/((1+x)*(1+t*x)*(t+z))
    assert sp.cancel(sp.residue(reduced/(1+z), t, -z)
                     -2/((1+x)*(1-z*x))) == 0
    # Blow up (z+1,t-1). On this chart the pole stays at v=-1.
    # This is the relative CD form, with z held fixed during t integration.
    pulled_back = sp.cancel(s*reduced.subs({z: s-1, t: 1+s*v}))
    expected = 2*s/((1+x)*(1+x*(1+s*v))*(1+v))
    assert sp.cancel(pulled_back-expected) == 0
    assert sp.cancel(sp.limit(pulled_back/s, s, 0)-2/((1+x)**2*(1+v))) == 0
    # A resolved meridian can approach the seam with nonzero winding in v.
    assert sp.cancel(sp.residue(pulled_back/s, v, -1)
                     -2/((1+x)*(1+x*(1-s)))) == 0


def lateral_detour(a, x, radius, side):
    """Integral of the bare dt kernel on a left-to-right semicircle."""
    def density(angle):
        offset = radius*exp(1j*angle)
        # dt=i*offset*dangle and the pole factor is offset.
        return 2j/((1+x)*(1+(a+offset)*x))

    return -side*quad(lambda angle: density(side*angle), 0, pi, complex_func=True,
                      epsabs=1e-11, epsrel=1e-11)[0]


def check_cd_cycle_seam():
    certify_moving_pole()
    for a in (0.5, 1., 2.):
        radius = a/4
        upper = a+radius*np.exp(1j*np.linspace(pi, 0, 65))
        lower = a+radius*np.exp(1j*np.linspace(-pi, 0, 65))
        loop = np.r_[upper, lower[-2::-1]]
        assert winding_number(loop, point=a) == -1
        for x in (0.3, 1., 4.):
            assert winding_number(loop, point=-1/x) == 0
            bare_jump = lateral_detour(a, x, radius, 1)-lateral_detour(a, x, radius, -1)
            np.testing.assert_allclose(bare_jump, -4j*pi/((1+x)*(1+a*x)),
                                       rtol=0, atol=1e-10)

        residue_integral = quad(lambda x: 1/((1+x)*(1+a*x)), 0, np.inf,
                                epsabs=1e-11, epsrel=1e-11)[0]
        expected = 1. if a == 1 else log(a)/(a-1)
        np.testing.assert_allclose(residue_integral, expected, atol=1e-11, rtol=0)
        stripped_jump = -4j*pi*residue_integral
        plus, minus = (integrated_projective(-a, log(a)+side*1j*pi) for side in (1, -1))
        np.testing.assert_allclose(plus-minus, (1-a)*stripped_jump, atol=1e-9, rtol=0)
        print(f"z=-{a:g}: detour winding=-1, scalar jump={plus-minus:.9f}, "
              f"jump without numerator factor={stripped_jump:.9f}")

    # At the seam the scalar jumps vanish. Dividing by z+1 keeps the two
    # germs distinct; evaluate the finite limit without subtracting pi^2.
    for side in (1, -1):
        errors = []
        for height in (1e-2, 1e-3, 1e-4):
            z = complex(-1, side*height)
            value = integrated_projective(z)/(1+z)
            errors.append(abs(value+side*2j*pi))
        assert errors[2] < errors[1] < errors[0] and errors[-1] < 0.001
        # The expanded benchmark derivative gives the same germ.
        np.testing.assert_allclose(2*complex_log(complex(-1, side*0.0))/(-1),
                                   -side*2j*pi, rtol=0, atol=1e-12)

    angles = np.linspace(0, -2*pi, 129)
    v_loop = -1+np.exp(1j*angles)/4
    assert winding_number(v_loop, point=-1) == -1
    for s in (0.1j, 0.001j):
        t_loop = 1+s*v_loop
        np.testing.assert_allclose(t_loop-1, s*v_loop, atol=1e-15, rtol=0)
        assert winding_number(t_loop, point=1-s) == -1

    print("PASS: CD parameter seam winding and the vanishing scalar jump")


if __name__ == "__main__":
    check_cd_cycle_seam()
