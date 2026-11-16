"""Compare the three cusp functionals on the Mellin family z**p.

The geometric capped-cycle period, the unit-circle mean and the leading
double pole of the dimensionally regulated frozen-loop integral are linear
in the finite weight F(z). This example evaluates all three on z**p and on
powers of log(z). See docs/GATE_6_FUNCTIONAL.md for the derivation and scope.
"""

from math import pi

import mpmath as mp
import numpy as np

from geocusp.benchmark import cusp_coefficient_from_mean

from .complex_cross_ratio import integrated_projective


def sinc(p):
    """Unit-circle mean of z**p on the principal sheet."""
    return mp.sinpi(p)/(mp.pi*p) if p else mp.mpf(1)


def circle_mean(weight, order=96):
    """Gauss-Legendre mean of weight(exp(i*phi)) over -pi < phi < pi."""
    nodes, weights = np.polynomial.legendre.leggauss(order)
    return weights @ np.array([weight(np.exp(1j*pi*t)) for t in nodes])/2


def lifted_cycle(alpha, phi, theta, radius):
    """Point (a, u, v) of the capped cycle on the weight-two cover."""
    a = radius*np.exp(1j*alpha)
    product = a*(np.exp(1j*phi)+1)
    u = np.sqrt(abs(product))*np.exp(1j*theta)
    return np.array([a, u, product/u])


def cycle_period(weight, radius=1.3, order=24, turns=6, step=1e-6):
    """Integrate weight(z)*2*da^du^dv/(a*(a-u*v)) over the capped cycle.

    The Jacobian is a central difference of the parameterization, so this
    does not assume the pulled-back form. The result carries the
    normalization (2*pi*i)**-3, with one turn of theta covering [0, pi].
    """
    nodes, weights = np.polynomial.legendre.leggauss(order)
    total = 0j
    for node, weight_phi in zip(nodes, weights):
        phi = pi*node
        for alpha in 2*pi*np.arange(turns)/turns:
            for theta in pi*(np.arange(turns)+.5)/turns:
                point = np.array([alpha, phi, theta])
                columns = []
                for k in range(3):
                    shift = np.zeros(3)
                    shift[k] = step
                    columns.append((lifted_cycle(*(point+shift), radius)
                                    - lifted_cycle(*(point-shift), radius))/(2*step))
                a, u, v = lifted_cycle(*point, radius)
                z = u*v/a-1
                form = 2/(a*(a-u*v))
                total += (weight_phi*pi*(2*pi/turns)*(pi/turns)*weight(z)*form
                          * np.linalg.det(np.array(columns).T))
    return total/(2j*pi)**3


def double_pole(p, eta, ratio=1):
    """eps**2*s*t times the frozen-loop integral of z**p at eps=-eta.

    Propagator powers are 1+p on legs 1,3 and 1-p on legs 2,4. One
    Mellin-Barnes integral remains after the Feynman parameters. The
    infrared regulator needs eps<0, and ratio=s/t with t=1.
    """
    eps = -mp.mpf(eta)
    center = -1-p+mp.mpf(eta)/2

    def integrand(y):
        w = center+1j*y
        return (mp.gamma(-w)*mp.gamma(2+eps+w)*mp.gamma(1+p+w)**2
                * mp.gamma(-1-p-eps-w)**2*mp.mpf(ratio)**(1+p+w))

    knots = [-mp.inf, -30, -5, -1, -10*eta, -eta, 0, eta, 10*eta, 1, 5, 30, mp.inf]
    integral = mp.quad(integrand, knots)/(2*mp.pi)
    prefactor = mp.rgamma(1+p)**2*mp.rgamma(1-p)**2*mp.rgamma(-2*eps)
    return (eps**2*prefactor*integral).real


def check_lorentzian_double_pole():
    mp.mp.dps = 25
    for p in (mp.mpf(0), mp.mpf("0.3"), mp.mpf("-0.6"), mp.mpf("0.85")):
        target = 4*sinc(p)
        errors = [abs(double_pole(p, eta)-target) for eta in (.1, .03, .01, .003)]
        assert all(new < old for old, new in zip(errors, errors[1:]))
        # The remainder is linear in the regulator. Remove it with one
        # Richardson step and compare with the closed form.
        fine, finer = double_pole(p, .004), double_pole(p, .002)
        extrapolated = 2*finer-fine
        assert abs(extrapolated-target) < 2e-4
        print(f"p={float(p):+.2f}: double pole -> {float(extrapolated):.6f}, "
              f"4*sin(pi*p)/(pi*p)={float(target):.6f}")

    # At p=0 this is the massless box. Its single pole is -2*log(s*t)/eps.
    slope = (double_pole(0, .004, 3)-double_pole(0, .004, 1))/(-.004)
    assert abs(slope+2*mp.log(3)) < .02
    print(f"p=0 single-pole slope in log(s/t): {float(slope):.5f}, box value={-2*np.log(3):.5f}")


def check_geometric_period():
    for p in (0., .3, -.6, .85, .4+.5j):
        exact = complex(sinc(mp.mpmathify(p)))
        np.testing.assert_allclose(circle_mean(lambda z: np.exp(p*np.log(z))), exact,
                                   rtol=0, atol=1e-12)
        for radius in (.7, 1.3):
            period = cycle_period(lambda z: np.exp(p*np.log(z)), radius)
            np.testing.assert_allclose(period, exact, rtol=0, atol=2e-8)
    print("Capped-cycle period of z**p equals sin(pi*p)/(pi*p) at both radii")

    # Powers of log(z) are derivatives in p at p=0.
    moments = {0: 1, 1: 0, 2: -pi**2/3, 3: 0, 4: pi**4/5}
    for power, moment in moments.items():
        period = cycle_period(lambda z: np.log(z)**power)
        np.testing.assert_allclose(period, moment, rtol=0, atol=3e-7)
        derivative = mp.diff(sinc, 0, power) if power else sinc(0)
        np.testing.assert_allclose(complex(derivative), moment, rtol=0, atol=1e-9)

    # The CD weight, evaluated by its own quadrature on the cycle.
    period = cycle_period(integrated_projective, order=32)
    np.testing.assert_allclose(period, 2*pi*pi/3, rtol=0, atol=3e-7)
    gamma = [cusp_coefficient_from_mean(1, cycle_period(lambda z: -1+0*z)),
             cusp_coefficient_from_mean(2, period)]
    np.testing.assert_allclose(gamma, [4, -4*pi*pi/3], rtol=0, atol=2e-6)
    print(f"Integrated CD period on the cycle={period.real:.9f}, "
          f"Gamma[1]={gamma[0].real:.6f}, Gamma[2]={gamma[1].real:.6f}")


def check_mellin_functional():
    check_geometric_period()
    check_lorentzian_double_pole()
    print("PASS: cycle period, circle mean and Lorentzian double pole agree on z**p")


if __name__ == "__main__":
    check_mellin_functional()
