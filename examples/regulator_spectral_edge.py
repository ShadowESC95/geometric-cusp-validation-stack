"""Small-h spectrum and the subleading small-regulator terms."""

from math import log, pi

import numpy as np
import sympy as sp

from .critical_regulator_kernel import critical_kernel
from .regulator_spectrum import spectral_density


def edge_coefficients(z):
    """Coefficients in sigma(h)=h*(A*log(1/h)+B)+O(h^2*log(1/h))."""
    scale = (1+z)**2/2
    endpoints = 1+1/z**2
    A = scale*endpoints
    B = scale*(endpoints*log(4/(1+z))+(2/z**2-1)*log(z)+6/z-endpoints)
    return A, B


def certify_edge_coefficients():
    z, q = sp.symbols("z q", positive=True)
    g = 1+(z-1)*q
    subtracted = (1/g**2-(1-q)-q/z**2)/(q*(1-q))
    primitive = -(z*z-1)*sp.log(g)/z**2+(z-1)/(z*g)
    assert sp.cancel(sp.diff(primitive, q)-subtracted) == 0
    endpoint_remainder = -sp.log(z)-1+2/z+sp.log(z)/z**2-1/z**2
    assert sp.simplify(primitive.subs(q, 1)-primitive.subs(q, 0)-endpoint_remainder) == 0

    # On the small-t branch q begins at (1+z)*h/4 and
    # 1-q ends at (1+z)*h/(4*z). The large-t branch adds 4/z.
    scale = (1+z)**2/2
    A = scale*(1+1/z**2)
    B = scale*((1+1/z**2)*sp.log(4/(1+z))+sp.log(z)/z**2+endpoint_remainder+4/z)
    for coefficient in (A, B):
        reciprocal = coefficient.subs(z, 1/z)
        assert sp.simplify(sp.expand_log(coefficient-reciprocal, force=True)) == 0
    np.testing.assert_allclose(
        edge_coefficients(1.5), [float(c.subs(z, sp.Rational(3, 2))) for c in (A, B)],
        rtol=0, atol=1e-12,
    )


def check_regulator_spectral_edge():
    certify_edge_coefficients()
    for z in (0.5, 1., 1.5, 5.):
        A, B = edge_coefficients(z)
        errors = []
        for h in (1e-4, 1e-5, 1e-6, 1e-7):
            residual = spectral_density(h, z)/h-A*log(1/h)
            errors.append(abs(residual-B))
        assert all(new < old for old, new in zip(errors, errors[1:]))
        assert errors[-1] < 7e-4
        print(f"z={z:g}: A={A:.9f}, B={B:.9f}, small-h estimate of B={residual:.9f}")

        boundary = critical_kernel(z, 0)
        real_errors, imaginary_remainders = [], []
        for eta in (1e-5, 1e-7, 1e-9):
            value = critical_kernel(z, eta)
            L = log(1/eta)
            real_coefficient = 2*(boundary-value.real)/(pi*eta)-A*L
            real_errors.append(abs(real_coefficient-B))
            imaginary_remainders.append(value.imag/eta-A*L**2/2-B*L)
        assert real_errors[-1] < min(real_errors[0], 0.001)
        # Subtracting both imaginary logarithms leaves a finite constant.
        assert abs(imaginary_remainders[-1]-imaginary_remainders[-2]) < 0.003
        print(f"z={z:g}, eta={eta:.0e}: real estimate of B={real_coefficient:.9f}, "
              f"imaginary remainder={imaginary_remainders[-1]:.9f}")
    print("PASS: low spectral edge and subleading regulator corrections")


if __name__ == "__main__":
    check_regulator_spectral_edge()
