"""Test a moving CD parameter square and the period of its secondary meridian."""

from cmath import exp, log
from math import pi

import numpy as np
import sympy as sp

from geocusp.contours import winding_number

from .cd_cycle_monodromy import meridian_period


def certify_compactified_form():
    x, y, X, Y, z = sp.symbols("x y X Y z")
    form = 2*(1+z)/((1+x)*(1+y)*(z*x+y))
    change = sp.Matrix([X/(1-X), Y/(1-Y)])
    determinant = change.jacobian([X, Y]).det()
    denominator = z*X*(1-Y)+Y*(1-X)
    pulled_back = form.subs(dict(zip((x, y), change)))*determinant
    assert sp.cancel(pulled_back-2*(1+z)/denominator) == 0

    # These opposite corners lie away from the moving divisor. Their
    # images under the same nonvanishing function have unequal winding.
    assert denominator.subs({X: 1, Y: 0}) == z
    assert denominator.subs({X: 0, Y: 1}) == 1

    # The residue along the moving divisor has a secondary pole where
    # that divisor runs to Y=infinity, which is y=-1.
    root = -z*X/(1-(1+z)*X)
    residue = sp.cancel(sp.residue(2*(1+z)/denominator, Y, root))
    assert sp.cancel(residue-2*(1+z)/(1-(1+z)*X)) == 0
    assert sp.cancel((X-1/(1+z))*residue) == -2
    U = sp.symbols("U")
    cleared = sp.cancel(U*denominator.subs(Y, 1/U))
    assert sp.cancel(cleared.subs(U, 0)-1+(1+z)*X) == 0
    at_infinity = sp.cancel(pulled_back.subs(Y, 1/U)*(-1/U**2))
    assert sp.cancel(U*at_infinity+2*(1+z)/cleared) == 0


def secondary_meridian_period(z, radius):
    """Integrate over positive X and clockwise Y meridians, in X,Y order."""
    pole = 1/(1+z)

    loop = pole+radius*np.exp(2j*pi*np.arange(256)/256)
    np.testing.assert_allclose(
        [winding_number(loop, point=p) for p in (pole, 0, 1)],
        [1, 0, 0], atol=1e-12, rtol=0,
    )
    # Evaluate the un-cancelled kernel too, on a representative inner
    # meridian. The coefficient of Y is nonzero for the entire X loop.
    def nested(theta, psi):
        X = pole+radius*exp(1j*theta)
        slope = 1-(1+z)*X
        Y0 = -z*X/slope
        radius_y = min(abs(Y0), abs(Y0-1))/4
        Y = Y0+radius_y*exp(-1j*psi)
        dX = 1j*(X-pole)
        dY = -1j*(Y-Y0)
        return 2*(1+z)*dX*dY/(z*X*(1-Y)+Y*(1-X))

    nodes, weights = np.polynomial.legendre.leggauss(12)
    angles = pi*(1+nodes)
    numeric = pi**2*(weights @ np.array([
        [nested(theta, psi) for psi in angles] for theta in angles
    ]) @ weights)
    return numeric


def check_cd_moving_contour():
    certify_compactified_form()
    circle = np.exp(2j*pi*np.arange(256)/256)
    assert winding_number(circle) == 1
    assert winding_number(np.ones_like(circle)) == 0
    print("Fixed corners of the moving square: divisor winding +1 versus 0")

    for angle in (-2.7, -pi/2, -0.2):
        z = exp(1j*angle)
        pole = 1/(1+z)
        radius = min(abs(pole), abs(pole-1))/4
        period = secondary_meridian_period(z, radius)
        np.testing.assert_allclose(period, -8*pi*pi, rtol=0, atol=1e-9)
        np.testing.assert_allclose(secondary_meridian_period(z, radius/2),
                                   period, rtol=0, atol=1e-9)
        continued = meridian_period(z)
        np.testing.assert_allclose(continued+period,
                                   4j*pi*(log(-z)+2j*pi), rtol=0, atol=1e-9)
        print(f"angle={angle:.3f}: secondary meridian period={period:.9f}")

    # The X winding is integer-valued on a pole-avoiding deformation.
    # A fixed winding can shift the branch by a constant but leaves
    # its change with arg(z) intact.
    z1, z2 = exp(-2j), exp(-0.5j)
    difference = meridian_period(z2)-meridian_period(z1)
    np.testing.assert_allclose(difference, -6*pi, rtol=0, atol=1e-9)
    print("PASS: moving-square winding obstruction and nested CD residue")


if __name__ == "__main__":
    check_cd_moving_contour()
