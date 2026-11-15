"""Transport the prescription and glue local Lorentzian contour patches."""

from itertools import combinations, product
from math import pi

import numpy as np
import sympy as sp
from numpy.polynomial.legendre import leggauss
from scipy.integrate import quad

from .geometry import (
    AB, EXTERNAL, EXTERNAL_PAIRS, H, INFINITY, bracket, lorentzian_line,
    mutual, numerator, regulated,
)


def check_transport():
    rho, eps = sp.symbols("rho eps", positive=True)
    a, d, A, D, u, v, U, V = sp.symbols("a d A D u v U V", real=True)
    G = sp.diag(1, 1/rho, 1, rho)
    assert G.det() == 1 and sp.simplify(G*H*G.T - H) == sp.zeros(4)

    def lift(line):
        return sp.diag(1, rho)*line*G

    ab, cd = lorentzian_line(a, d, u+sp.I*v), lorentzian_line(A, D, U+sp.I*V)
    ab_scaled, cd_scaled = lift(ab), lift(cd)
    powers = (2, 2, 0, 0)
    for (i, j), power in zip(EXTERNAL_PAIRS, powers):
        for original, moved in ((ab, ab_scaled), (cd, cd_scaled)):
            assert sp.expand(bracket(moved, i, j) - rho**power*bracket(original, i, j)) == 0
    assert sp.expand(mutual(ab_scaled, cd_scaled) - rho**2*mutual(ab, cd)) == 0
    assert sp.expand(numerator(ab_scaled, cd_scaled) - rho**2*numerator(ab, cd)) == 0
    coords = sp.Matrix([A, rho**2*D, rho*U, rho*V])
    assert coords.jacobian([A, D, U, V]).det() == rho**4

    cd = lorentzian_line(-2, sp.Rational(4, 5), -1)
    moved_cd = lift(cd)
    for i, (line, moved, power) in enumerate(zip(
        EXTERNAL+[AB], EXTERNAL+[lift(AB)], powers+(2,),
    )):
        assert sp.cancel(regulated(moved_cd, moved, INFINITY*G, eps)/rho**power
                         - regulated(cd, line, INFINITY, eps)) == 0
        fixed_shift = -sp.im(regulated(moved_cd, moved, INFINITY, eps)/rho**power)/eps
        limit = sp.limit(fixed_shift, rho, 0, dir="+")
        assert (limit == sp.oo) == (i in (0, 1, 4))
    print("All five i0 prescriptions transport exactly. Three fixed-frame shifts diverge.")

    s = sp.symbols("s", positive=True)

    def projective_limit(line, index):
        minors = [line[:, [i,j]].det() for i,j in combinations(range(4), 2)]
        return [sp.limit(value/minors[index], rho, 0, dir="+") for value in minors]

    fixed = lorentzian_line(a, s, 0)
    edge = lorentzian_line(a, s/rho**2, 0)
    assert projective_limit(lift(fixed), 0) == [1, 0, 0, -a, 0, 0]
    assert projective_limit(edge, 2) == [0, 0, 1, 0, 0, a]
    assert lift(edge) == fixed
    print("A chart-edge family stays away from the cusp after transport.")


def check_edge_residue():
    A, T, u, v, eps = sp.symbols("A T u v eps", real=True)
    b = u + sp.I*v
    cd = sp.Matrix([[1, b, A, 0], [0, T, 2*sp.conjugate(b), 1]])
    assert sp.simplify(cd*H*cd.conjugate().T) == sp.zeros(2)
    u12 = cd[:, :2].inv()*cd
    coords = sp.Matrix([u12[0, 2], u12[1, 3], sp.re(u12[0, 3]), sp.im(u12[0, 3])])
    jac = sp.factor(coords.jacobian([A, T, u, v]).det())

    def form(line):
        factors = [bracket(line, i, j) for i,j in EXTERNAL_PAIRS]
        return numerator(AB, line)/(mutual(AB, line)*sp.prod(factors))

    assert jac == -1/T**4 and sp.cancel(form(u12)*jac + form(cd)) == 0
    point = {A: -1, u: sp.Rational(1, 3), v: 0}
    residue = sp.residue(sp.cancel((form(u12)*jac).subs(point)), T, 0)
    shift = regulated(cd, EXTERNAL[2], INFINITY, eps).subs({**point, T: 0})
    assert residue == sp.Rational(27, 5) and shift == 8*sp.I*eps/9
    print(f"Oriented edge residue: {residue}. Regulated denominator: {shift}.")


def chart_density(coordinate, A, u, v, eps, *, reciprocal=False):
    """One oriented T interval, or its d=1/T representation at fixed (A,u,v)."""
    scale = coordinate if reciprocal else 1
    t = 1/coordinate if reciprocal else coordinate
    b_squared = 2*(u*u+v*v)
    edge = -A*t + b_squared
    factors = scale * np.array([A, 1, t, edge, edge+A+t-4*u+1])
    infinity = scale*(edge+A+t-1)
    den = factors - 0.5j*eps*infinity*np.array([1, 1, 1, -1, 2])
    numerator = -4*u*scale
    jac = coordinate**2 if reciprocal else 1
    return -jac*numerator/np.prod(den)


def normal_integrals(A, u, v, eps):
    def integrate(lo, hi, reciprocal=False):
        return quad(lambda T: chart_density(T, A, u, v, eps, reciprocal=reciprocal),
                    lo, hi, complex_func=True, epsabs=1e-9, epsrel=1e-9,
                    limit=250, points=[0] if lo < 0 < hi else None)[0]

    width = 0.05
    whole = integrate(-width, width)
    splits = []
    for cut in (0.01, 0.02):
        splits.append(integrate(-1/width, -1/cut, reciprocal=True)
                      + integrate(-cut, cut)
                      + integrate(1/cut, 1/width, reciprocal=True))
    return np.array([whole, *splits])


def check_chart_matching():
    for eps in (1e-3, 1e-4, 1e-5):
        whole, *splits = normal_integrals(-1, 1/3, 0, eps)
        np.testing.assert_allclose(splits, whole, rtol=0, atol=1e-8)
        print(f"Fixed slice eps={eps:.0e}: integral={whole:.9f}")
    assert abs(whole.imag + 27*pi/5) < 2e-3

    errors = []
    for order, eps in ((3, 1e-3), (3, 1e-4), (3, 1e-5), (5, 1e-5)):
        nodes, weights = leggauss(order)
        axes = [zip((lo+hi)/2+(hi-lo)/2*nodes, (hi-lo)/2*weights)
                for lo, hi in ((-1.05,-0.95), (0.30,0.36), (-0.03,0.03))]
        integral = np.zeros(3, dtype=complex)
        target = 0.
        for (A, wa), (u, wu), (v, wv) in product(*axes):
            weight = wa*wu*wv
            b_squared = 2*(u*u+v*v)
            mutual_at_edge = b_squared+A-4*u+1
            assert b_squared+A-1 < 0 and mutual_at_edge < 0
            target -= weight*pi*4*u/(A*b_squared*mutual_at_edge)
            integral += weight*normal_integrals(A, u, v, eps)
        np.testing.assert_allclose(integral[1:], integral[0], rtol=0, atol=1e-9)
        errors.append(abs(integral[0].imag - target))
        print(f"Patch n={order} eps={eps:.0e}: Im={integral[0].imag:.12f}, "
              f"residue target={target:.12f}")
    assert errors[2] < errors[1] < errors[0] and errors[-1] < 1e-6


if __name__ == "__main__":
    check_transport()
    check_edge_residue()
    check_chart_matching()
