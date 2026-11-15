"""Exact geometry checks and the local pole diagnostics on two slices."""

from cmath import log, sqrt
from math import acos, cos, pi

import numpy as np
import sympy as sp
from scipy.integrate import quad

from geocusp.contours import winding_number
from geocusp.linking import winding
from geocusp.negative_geometry import (
    iterated_residue_order_scan,
    symbolic_two_loop_integrand,
    two_loop_blowup_leading_coefficient,
)
from geocusp.residue import (
    blowup_intersection_yw,
    strict_cusp_blowup_zero_residue,
    strict_cusp_weighted_blowup,
)
from geocusp.twistor import all_four_brackets, line_from_xyzw

Z = sp.eye(4)
J = sp.diag(1, -2)
H = sp.BlockMatrix([[sp.zeros(2), sp.I*J], [-sp.I*J, sp.zeros(2)]]).as_explicit()
INFINITY = Z[:2, :] + Z[2:, :]
EXTERNAL = [Z[[i, j], :] for i, j in ((0, 1), (1, 2), (2, 3), (0, 3))]
EXTERNAL_PAIRS = ((1, 2), (2, 3), (3, 4), (1, 4))
AB = sp.Matrix([[1, 0, -1, -1], [0, 1, 2, 1]])


def mutual(L, R):
    return L.col_join(R).det()


def bracket(L, i, j):
    return mutual(L, Z[[i-1, j-1], :])


def numerator(L, R):
    return -(bracket(L, 1, 3)*bracket(R, 2, 4)
             + bracket(L, 2, 4)*bracket(R, 1, 3))


def regulated(L, R, infinity, eps):
    return mutual(L, R) - sp.I*eps*mutual(L, infinity)*mutual(R, infinity)/2


def check_local_geometry():
    x, z0, rho, eta = sp.symbols("x z_0 rho eta", nonzero=True)
    residue, ratio = blowup_intersection_yw()
    assert residue == 1/(x*eta*z0) and ratio == eta*z0/x
    assert strict_cusp_blowup_zero_residue()[1] == 0
    weighted, _ = strict_cusp_weighted_blowup()
    A, z, D = sp.symbols("a z D", nonzero=True)
    assert sp.cancel(weighted + 1/(A*z*D)) == 0
    C = -A*D*(1+z)
    weighted_line = sp.Matrix([[1, 0, A, -rho], [0, 1, rho*C, rho**2*D]])
    assert sp.cancel(bracket(weighted_line, 1, 2)*bracket(weighted_line, 3, 4)
                     / (bracket(weighted_line, 2, 3)*bracket(weighted_line, 1, 4)) - z) == 0

    y, w = rho*z*x/z0, rho
    jac = sp.Matrix([x, y, z0, w]).jacobian([rho, x, z, z0]).det()
    sign = sp.cancel(rho*x*z*z0*jac/(x*y*z0*w))
    assert sign == -1
    circle = np.exp(2j*np.pi*np.arange(1024)/1024)
    for measure in (winding_number, winding):
        np.testing.assert_allclose(
            [measure(circle), measure(circle[::-1]), measure(circle, point=2)],
            [1, -1, 0], rtol=0, atol=1e-12,
        )
    period = float(sign) * np.prod([winding_number(r*circle) for r in (1.3, 1., 0.9)])
    np.testing.assert_allclose(period, -1, rtol=0, atol=1e-12)
    print(f"Radial-first signed product period: {period:g}")

    np.testing.assert_allclose(
        list(all_four_brackets(line_from_xyzw(2., 3., 4., 5.)).values()),
        [3, 5, 4, 2, -22, -1],
    )
    variables, form, X12, n12 = symbolic_two_loop_integrand()
    x1, y1, z1, w1, x2, y2, z2, w2 = variables
    L = sp.Matrix([[1, x1, 0, -y1], [0, z1, 1, w1]])
    R = sp.Matrix([[1, x2, 0, -y2], [0, z2, 1, w2]])
    assert sp.expand(mutual(L, R) - X12) == 0
    assert sp.expand(numerator(L, R) + n12) == 0
    rescaled, lead = two_loop_blowup_leading_coefficient()
    radial = sp.Symbol("rho")
    assert sp.cancel(radial**4 * rescaled - lead) == 0
    scan = iterated_residue_order_scan()
    assert scan["n_orders"] == 24 and scan["n_unique_results"] == 4
    print("Two-loop homogeneity: -4. Residue orders: 24 with four distinct values.")

    fixed = sp.cancel((x1*y1*z1*w1*form).subs(
        {x1: 1, z1: 1, x2: 2, z2: 3, y1: radial*z, w1: radial}))
    A, B = sp.symbols("A B", positive=True)
    assert not sp.limit(fixed, radial, 0).has(z)
    inner = sp.cancel(radial**2*fixed.subs({y2: radial*A, w2: radial*B}))
    assert not inner.has(radial) and inner.has(z)
    print("Fixed CD limit loses z. The scaled CD kernel retains it.")


def check_complex_ray():
    rho, r = sp.symbols("rho r", positive=True)
    z = sp.symbols("z")

    def line(scale, a, D, cross):
        C = -a*D*(1 + cross)
        return sp.Matrix([[1, 0, a, -scale], [0, 1, scale*C, scale**2*D]])

    L, R = line(rho, -1, 1, z), line(rho*r, -2, 2, 1)
    Q = sp.cancel(mutual(L, R)/rho**2)
    N = sp.cancel(numerator(L, R)/rho**2)
    assert sp.expand(Q - (6*r*r - (z+9)*r + z+2)) == 0
    assert sp.expand(N + r*(z+9)) == 0
    roots = sp.solve(Q.subs(z, 1), r)
    assert len(roots) == 2 and all(root.is_positive for root in roots)
    assert all(sp.simplify(N.subs({z: 1, r: root})) != 0 for root in roots)
    # The external null square leaves two complex entries in its Hermitian form.
    p_real, p_imag, q_real, q_imag = sp.symbols("p_r p_i q_r q_i", real=True)
    J_general = sp.diag(p_real+sp.I*p_imag, q_real+sp.I*q_imag)
    H_general = sp.BlockMatrix(
        [[sp.zeros(2), J_general], [J_general.conjugate().T, sp.zeros(2)]],
    ).as_explicit()
    conditions = []
    for candidate in (L.subs(z, 1), R):
        for entry in candidate*H_general*candidate.conjugate().T:
            conditions.extend(entry.as_real_imag())
    assert sp.linsolve(conditions, (p_real, p_imag, q_real, q_imag)) == sp.FiniteSet((0, 0, 0, 0))

    angles = np.linspace(-pi, pi, 2049)
    samples = np.exp(1j*angles)
    discriminant = np.sqrt((samples - 3)**2 + 24)
    pair = np.array([(samples + 9 - discriminant)/12,
                     (samples + 9 + discriminant)/12])
    off_axis = np.abs(np.sin(angles)) > 1e-10
    assert np.all(pair[:, off_axis].imag*np.sin(angles[off_axis]) > 0)
    np.testing.assert_allclose(pair[:, 0], pair[:, -1], rtol=0, atol=1e-12)
    assert min(abs(pair[1] - pair[0])) > 0.8

    def ray_integral(cross, shift):
        small, large = sorted(
            np.roots([6, -(cross+9), cross+2+shift]), key=lambda root: root.real)
        residue = -(cross+9)/(6*(small-large))
        return residue*(log(-large) - log(-small))

    eps = 1e-7
    upward = sorted(np.roots([6, -10-1j*eps, 3+1j*eps]), key=lambda root: root.real)
    mixed = sorted(np.roots([6, -10, 3+1j*eps]), key=lambda root: root.real)
    np.testing.assert_array_equal(np.sign(np.imag(upward)), [1, 1])
    np.testing.assert_array_equal(np.sign(np.imag(mixed)), [1, -1])
    jump = ray_integral(1, 1j*eps) - ray_integral(1+1j*eps, 0j)
    np.testing.assert_allclose(jump, 10j*pi/np.sqrt(7), rtol=0, atol=2e-6)
    print(f"Complex-ray poles: {sorted(map(float, roots))}. Bypass jump: {jump:.9f}")


def lorentzian_line(a, d, b):
    return sp.Matrix([[1, 0, a, b], [0, 1, -2*sp.conjugate(b), d]])


def affine(line):
    U, V = line[:, :2], line[:, 2:]
    scale = sp.factor((V-U).det())
    X = sp.simplify((V-U).inv()*U*J)
    assert X.is_hermitian
    return scale, X


def check_lorentzian_poles():
    rho, r, theta = sp.symbols("rho r theta", positive=True)
    L = lorentzian_line(-1, rho**2, -rho)
    control = lorentzian_line(-2, rho**2*r**2/2, -rho*r)
    R = lorentzian_line(-2, 4*rho**2*r**2/5, -rho*r*sp.exp(sp.I*theta))
    for line in (L, R, control):
        assert sp.simplify(line*H*line.conjugate().T) == sp.zeros(2)
    q_control = sp.cancel(mutual(L, control)/rho**2)
    assert q_control.subs(r, sp.Rational(4, 3)) == sp.Rational(1, 3)
    Q = sp.simplify(sp.expand_complex(mutual(L, R)/rho**2))
    assert sp.simplify(Q - (6*r*r/5 - 4*r*sp.cos(theta) + 3)) == 0

    def poles(phase, cross, eps=0.):
        return np.roots([2 - 1/(1+cross), -4*cos(phase), 3-0.5j*eps])

    collision = acos(3/np.sqrt(10))
    np.testing.assert_allclose(np.sort_complex(poles(0, 0.25)),
                               [(10-np.sqrt(10))/6, (10+np.sqrt(10))/6])
    assert np.all(abs(poles(collision + 0.05, 0.25).imag) > 0.1)
    shifted = sorted(poles(0, 0.25, 1e-7), key=lambda value: value.real)
    np.testing.assert_array_equal(np.sign(np.imag(shifted)), [-1, 1])
    print(f"Lorentzian pole collision at theta={collision:.9f}")

    # At this r, the mutual pole is inside the phase integral.
    radial = 5/3
    p, amp = 3 + 6*radial**2/5, 4*radial
    for eps in (1e-2, 1e-4):
        q = p - 0.5j*eps
        expected = 1 - q/(sqrt(q-amp)*sqrt(q+amp))
        actual = quad(lambda phase: -amp*cos(phase)/(q-amp*cos(phase)),
                      0, pi, points=[acos(p/amp)], complex_func=True,
                      epsabs=2e-9, epsrel=2e-9, limit=250)[0]/pi
        np.testing.assert_allclose(actual, expected, rtol=0, atol=2e-8)
    np.testing.assert_allclose(actual.imag, -p/np.sqrt(amp**2-p**2),
                               rtol=0, atol=1e-6)

    A, z = sp.symbols("A z", positive=True)
    d = 2*r*r/(A*(1+z))
    X = sp.Matrix([(-A-2*rho**2*d)/2, 2*rho*r*sp.cos(theta),
                   -2*rho*r*sp.sin(theta), (-A+2*rho**2*d)/2])
    jac = sp.trigsimp(X.jacobian([A, r, z, theta]).det())
    assert sp.cancel(jac/(rho**4*A*d*(2*r*r-A*d)) - 2/(A*r*z)) == 0

    for z, expected in ((0.25, -4*pi), (1., 0.)):
        k = 2 - 1/(1+z)
        discriminant = 16 - 12*k
        value = 0.
        if discriminant > 0:
            lo, hi = [((4+sign*np.sqrt(discriminant))/(2*k))**2 for sign in (-1, 1)]
            value = quad(lambda phi: -(3/((lo+hi+(hi-lo)*cos(phi))/2)+k)/(2*z*k),
                         0, pi, epsabs=1e-10, epsrel=1e-10)[0]
        np.testing.assert_allclose(value, expected, rtol=0, atol=1e-9)
        print(f"z_CD={z:g}: partial mutual-pole imaginary integral={value:.12f}")

    R = R.subs(theta, 0)
    scale_cd, X_cd = affine(R)
    signs = []
    lines = EXTERNAL + [L]
    for line in lines:
        scale, X = affine(line)
        assert sp.cancel((X_cd-X).det() + 2*mutual(R, line)/(scale_cd*scale)) == 0
        signs.append(sp.sign((-scale_cd*scale/2).subs({rho: sp.Rational(1, 2), r: 1})))
    assert signs == [-1, -1, -1, 1, -1]
    points = [affine(line.subs(rho, sp.Rational(1, 2)))[1] for line in lines]
    alpha = [sp.Rational(n, 100) for n in (19, 21, 19, 21, 20)]
    delta = -sum(alpha[i]*alpha[j]*(points[i]-points[j]).det()
                 for i in range(5) for j in range(i+1, 5))
    assert delta == 0 and sum(alpha) == 1 and all(a > 0 for a in alpha)
    print(f"Finite-frame propagator signs: {signs}. Interior parameter zero verified.")


if __name__ == "__main__":
    check_local_geometry()
    check_complex_ray()
    check_lorentzian_poles()
