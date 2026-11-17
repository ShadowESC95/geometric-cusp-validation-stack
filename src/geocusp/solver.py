"""Local cusp functional solver: numerics and command line.

Four linear functionals turn a finite weight F(z) into a cusp coefficient:
the unit-circle mean, the period over the capped cycle on the weighted
divisor, the leading double pole of the dimensionally regulated frozen-loop
integral, and the Lorentzian edge chain (the coefficient of the collinear
singularity at a light-like edge, divided by -2*pi*i). See
docs/GATE_6_FUNCTIONAL.md, docs/GATE_6_CLOSURE.md and docs/GATE_6_BOUNDS.md.

Run ``python -m geocusp.solver`` to open the desktop window, or use the
``functional``, ``edge``, ``two-loop``, ``contour``, ``cd`` and ``rim``
subcommands in a terminal. Everything runs in this process: there is no
server and no network access.
"""
from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from math import factorial, pi

import numpy as np
from scipy.integrate import simpson
from scipy.special import loggamma, rgamma

from .benchmark import cusp_coefficient_from_mean, weak_coefficients

ROTATION = np.exp(1j*pi/4)
CAP_RADIUS = 0.06
MAX_LOG_POWER = 4


# ---------------------------------------------------------------- weights

@dataclass(frozen=True)
class Weight:
    """F(z) = sum_k logs[k]*log(z)**k + amplitude*z**exponent on the principal sheet."""

    logs: tuple[float, ...] = (0.0,)
    amplitude: float = 0.0
    exponent: float = 0.0

    def __post_init__(self):
        if len(self.logs) > MAX_LOG_POWER+1:
            raise ValueError(f"at most log(z)**{MAX_LOG_POWER} is supported")
        if not -1 < self.exponent < 1:
            raise ValueError("the exponent p must satisfy -1 < p < 1")

    def __call__(self, z):
        log_z = np.log(np.asarray(z, dtype=complex))
        value = sum(c*log_z**k for k, c in enumerate(self.logs))
        return value+self.amplitude*np.exp(self.exponent*log_z)

    def of_log(self, L):
        """The weight as a function of L = log z, for mpmath arguments on any sheet."""
        import mpmath as mp

        value = sum(mp.mpf(c)*L**k for k, c in enumerate(self.logs) if c)
        if self.amplitude:
            value += mp.mpf(self.amplitude)*mp.exp(mp.mpf(self.exponent)*L)
        return value

    def describe(self) -> str:
        names = ["", "*log(z)", "*log(z)**2", "*log(z)**3", "*log(z)**4"]
        terms = [f"{c:g}{names[k]}" for k, c in enumerate(self.logs) if c]
        if self.amplitude:
            terms.append(f"{self.amplitude:g}*z**{self.exponent:g}")
        return " + ".join(terms).replace("+ -", "- ") or "0"


def sinc(p):
    """Unit-circle mean of z**p, sin(pi*p)/(pi*p)."""
    return float(np.sinc(p))


def log_moment(power: int) -> float:
    """Unit-circle mean of log(z)**power: zero for odd powers."""
    if power % 2:
        return 0.0
    return (-1)**(power//2)*pi**power/(power+1)


def closed_mean(weight: Weight) -> float:
    """Exact circle mean of the weight."""
    logs = sum(c*log_moment(k) for k, c in enumerate(weight.logs))
    return logs+weight.amplitude*sinc(weight.exponent)


# ---------------------------------------------------------------- functionals

def circle_mean(weight, order: int = 48) -> complex:
    """Gauss-Legendre mean of weight(exp(i*phi)) over -pi < phi < pi."""
    nodes, weights = np.polynomial.legendre.leggauss(order)
    return complex(weights @ weight(np.exp(1j*pi*nodes))/2)


def lifted_cycle(alpha, phi, theta, radius):
    """Point (a, u, v) of the capped cycle on the weight-two cover."""
    a = radius*np.exp(1j*alpha)
    product = a*(np.exp(1j*phi)+1)
    u = np.sqrt(abs(product))*np.exp(1j*theta)
    return a, u, product/u


def cycle_period(weight, radius: float = 1.3, order: int = 32, turns: int = 6,
                 step: float = 1e-6) -> complex:
    """Integrate weight(z)*2*da^du^dv/(a*(a-u*v)) over the capped cycle.

    The Jacobian is a central difference of the parameterization. The result
    carries the normalization (2*pi*i)**-3, with theta covering [0, pi].
    """
    if radius <= 0:
        raise ValueError("the cycle radius must be positive")
    nodes, weights = np.polynomial.legendre.leggauss(order)
    alpha, phi, theta = np.meshgrid(2*pi*np.arange(turns)/turns, pi*nodes,
                                    pi*(np.arange(turns)+.5)/turns, indexing="ij")
    point = [alpha, phi, theta]
    columns = []
    for k in range(3):
        plus = [value+step*(k == j) for j, value in enumerate(point)]
        minus = [value-step*(k == j) for j, value in enumerate(point)]
        columns.append([(p-m)/(2*step) for p, m in zip(lifted_cycle(*plus, radius),
                                                         lifted_cycle(*minus, radius))])
    (a0, u0, v0), (a1, u1, v1), (a2, u2, v2) = columns
    determinant = (a0*(u1*v2-u2*v1)-a1*(u0*v2-u2*v0)+a2*(u0*v1-u1*v0))
    a, u, v = lifted_cycle(alpha, phi, theta, radius)
    form = 2/(a*(a-u*v))
    measure = weights[None, :, None]*pi*(2*pi/turns)*(pi/turns)
    return complex(np.sum(measure*weight(u*v/a-1)*form*determinant)/(2j*pi)**3)


def double_pole(p: float, eta: float, ratio: float = 1.0, points: int = 1400) -> float:
    """eps**2*s*t times the frozen-loop integral of z**p at eps=-eta.

    Legs 1,3 carry power 1+p and legs 2,4 carry 1-p. One Mellin-Barnes
    integral remains. The infrared regulator needs eps<0; ratio=s/t, t=1.
    As eta tends to zero the value tends to 4*sin(pi*p)/(pi*p).
    """
    if not 0 < eta <= .5:
        raise ValueError("the regulator must satisfy 0 < eta <= 0.5")
    if not -1 < p < 1:
        raise ValueError("the exponent p must satisfy -1 < p < 1")
    eps = -eta
    extent = np.arcsinh(30/eta)
    t = np.linspace(-extent, extent, points+1)
    w = -1-p+eta/2+1j*eta*np.sinh(t)
    log_integrand = (loggamma(-w)+loggamma(2+eps+w)+2*loggamma(1+p+w)
                     + 2*loggamma(-1-p-eps-w)+(1+p+w)*np.log(ratio))
    integral = simpson((np.exp(log_integrand)*eta*np.cosh(t)).real, x=t)/(2*pi)
    return float(eps**2*rgamma(1+p)**2*rgamma(1-p)**2*rgamma(-2*eps)*integral)


def lorentzian(weight: Weight, eta: float) -> float:
    """Regulated Lorentzian functional of the weight, normalized by 4.

    Powers of log(z) are derivatives in p at p=0. The regulated pole is even
    in p at s=t, so an even polynomial fit supplies those derivatives.
    """
    total = 0.0
    if any(weight.logs):
        nodes = np.linspace(0, .6, 7)
        values = np.array([double_pole(p, eta) for p in nodes])/4
        even = np.linalg.solve(np.vander(nodes**2, increasing=True), values)
        for power, coefficient in enumerate(weight.logs):
            if coefficient and power % 2 == 0:
                total += coefficient*factorial(power)*even[power//2]
    if weight.amplitude:
        total += weight.amplitude*double_pole(weight.exponent, eta)/4
    return float(total)


REGIONS = {"--": (-1.0, -1.7), "+-": (1.3, -1.0), "-+": (-1.0, 0.8)}
TWO_LOOP_MEAN = -11*pi**4/15


def edge_chain(weight, cusp: int = 2, edge: str = "next", position: float = .4,
               region: str = "--", ratio: float = 1e-9, digits: int = 15) -> dict:
    """Lorentzian edge chain: the collinear singularity at a point of a null line.

    ``weight`` is a Weight or the string "two-loop". ``position`` runs along
    the null line through the cusp and its neighbour; the edge segment is
    0 < position < 1. Returns the jump of the transverse kernel across the
    collinear invariant divided by -2*pi*i, which equals the unit-circle mean
    of the weight on the segment and zero off it.
    """
    import mpmath as mp

    from . import lorentzian_edges as le

    if cusp not in (1, 2, 3, 4):
        raise ValueError("cusp must be 1, 2, 3 or 4")
    if edge not in ("next", "prev"):
        raise ValueError("edge must be 'next' or 'prev'")
    if region not in REGIONS:
        raise ValueError("region must be one of '--', '+-', '-+' (the signs of s and t)")
    if position in (0, 1):
        raise ValueError("position 0 and 1 are the cusps; choose a point of the line off them")
    two_loop = isinstance(weight, str)
    if two_loop and weight != "two-loop":
        raise ValueError("the only named weight is 'two-loop'")
    f = le.f2_cover if two_loop else weight.of_log
    geometry = le.Cusp(cusp, *REGIONS[region])
    with mp.workdps(digits):
        jump = geometry.edge_jump(f, edge, position, mp.mpf(ratio))
        constant = geometry.edge_jump(lambda L: 1, edge, position, mp.mpf(ratio))
        value = complex(jump/(-2j*mp.pi))
        unit = complex(constant/(-2j*mp.pi))
        sheets = sorted({k for a in (1, -1) for b in (1, -1)
                         for k in geometry.sheet_sequence(a*1e-3*abs(geometry.D),
                                                          b*1e-7*abs(geometry.D))})
    on_segment = 0 < position < 1
    mean = TWO_LOOP_MEAN if two_loop else closed_mean(weight)
    expected = mean if on_segment else 0.
    return {
        "weight": "two-loop weight F2" if two_loop else weight.describe(),
        "cusp": cusp, "edge": edge, "position": position, "region": region,
        "on_segment": on_segment,
        "edge_chain": value.real, "edge_chain_imag": value.imag,
        "constant_weight": unit.real,
        "circle_mean": mean,
        "expected": expected,
        "agree": bool(abs(value-expected) < 1e-5*max(1., abs(mean))),
        "sheets": sheets,
    }


def edge_scan(weight, positions, cusp: int = 2, edge: str = "next", region: str = "--") -> list:
    """Edge chain at several positions along one null line."""
    return [edge_chain(weight, cusp, edge, x, region, ratio=1e-8, digits=13) for x in positions]


def solve_functional(weight: Weight, loop_order: int = 2, radius: float = 1.3,
                     eta: float = .01, curve: bool = False, edge: bool = False) -> dict:
    """Evaluate the functionals and the resulting cusp coefficient.

    With ``edge=True`` the Lorentzian edge chain is evaluated as well (about a
    second), at the midpoint region of an edge at cusp 2.
    """
    if loop_order < 1:
        raise ValueError("loop_order must be positive")
    exact = closed_mean(weight)
    circle = circle_mean(weight)
    cycle = cycle_period(weight, radius)
    regulated, halved = lorentzian(weight, eta), lorentzian(weight, eta/2)
    extrapolated = 2*halved-regulated
    scale = max(1., abs(exact))
    reference = weak_coefficients().get(loop_order)
    gamma = float(np.real(cusp_coefficient_from_mean(loop_order, exact)))
    agree = abs(cycle-exact) < 1e-5*scale and abs(extrapolated-exact) < 2e-3*scale
    result = {
        "weight": weight.describe(),
        "closed_mean": exact,
        "circle_mean": circle.real,
        "cycle_period": cycle.real,
        "cycle_period_imag": cycle.imag,
        "lorentzian": regulated,
        "lorentzian_extrapolated": extrapolated,
        "loop_order": loop_order,
        "gamma": gamma,
        "reference": reference,
        "matches_reference": bool(reference is not None and abs(gamma-reference) < 5e-3),
        "radius": radius,
        "eta": eta,
    }
    if edge:
        chain = edge_chain(weight)
        result["edge_chain"] = chain["edge_chain"]
        agree = agree and chain["agree"]
    result["agree"] = bool(agree)
    if curve:
        etas = (.1, .06, .03, .02, .01, .006, .003)
        result["curve"] = [{"eta": e, "value": lorentzian(weight, e)} for e in etas]
    return result


def solve_two_loop(edge: bool = True, position: float = .4) -> dict:
    """The published two-loop weight F2 and the three-loop cusp coefficient.

    F2 is not a polynomial in log z, so the Mellin-Barnes and capped-cycle
    routines for the Weight family do not apply; the circle mean is computed
    directly and, with ``edge=True``, so is the Lorentzian edge chain on the
    continued sheets (about a minute).
    """
    import mpmath as mp

    from . import lorentzian_edges as le

    with mp.workdps(18):
        circle = complex(le.circle_mean(le.f2_cover))
        at_minus_one = complex(le.f2_regular_part(0))
    reference = weak_coefficients()[3]
    gamma = float(np.real(cusp_coefficient_from_mean(3, circle)))
    result = {
        "weight": "two-loop weight F2 (arXiv:2112.06956, Eq. 2.23)",
        "closed_mean": TWO_LOOP_MEAN,
        "circle_mean": circle.real,
        "value_at_z_minus_one": at_minus_one.real,
        "regular_at_z_minus_one": True,
        "loop_order": 3,
        "gamma": gamma,
        "reference": reference,
        "matches_reference": bool(abs(gamma-reference) < 1e-6),
    }
    agree = abs(circle-TWO_LOOP_MEAN) < 1e-9
    if edge:
        chain = edge_chain("two-loop", position=position)
        off = edge_chain("two-loop", position=-position)
        result["edge_chain"] = chain["edge_chain"]
        result["edge_chain_off_segment"] = off["edge_chain"]
        agree = agree and chain["agree"] and off["agree"]
    result["agree"] = bool(agree)
    return result


def contour_survey(points: int = 5000, seed: int = 0) -> dict:
    """Survey the explicit integration contour of docs/GATE_6_CONTOUR.md."""
    from . import feynman_contour as fc

    if points < 100:
        raise ValueError("use at least 100 sample points")
    contour = fc.Contour()
    result = {key: float(value) for key, value in fc.survey(contour, points, seed).items()}
    result["near_count"] = int(result["near_count"])
    result["points"] = points
    result["feynman_side_of_every_cone"] = bool(result["max_im_near"] <= 0)
    result["zeros_only_on_edges"] = bool(result["min_ratio"] > 0)
    return result


# ---------------------------------------------------------------- CD weight

def cd_weight(z: complex, points: int = 3000) -> complex:
    """Reduced projective CD integral at complex z, off the branch point z=0."""
    z = complex(z)
    if abs(z) < 1e-9:
        raise ValueError("z=0 is a branch point")
    t = np.linspace(-45, 45, points+1)
    y = np.exp(t)
    offset = (y-z)/z
    small, middle = abs(offset) < 1e-3, abs(offset) < .5
    safe = np.where(small, 1, y-z)
    near = np.log1p(np.where(middle, offset, 0))
    quotient = np.where(middle, near, t-np.log(z))/safe
    quotient = np.where(small, (1-offset/2+offset**2/3)/z, quotient)
    return complex(2*(1+z)*simpson(quotient*y/(1+y), x=t))


def solve_cd(z: complex) -> dict:
    value = cd_weight(z)
    log_z = np.log(complex(z))
    expected = log_z**2+pi**2
    return {"z": [z.real, z.imag], "value": [value.real, value.imag],
            "closed_form": [expected.real, expected.imag],
            "error": abs(value-expected), "log_z": [log_z.real, log_z.imag],
            "modulus": abs(z), "on_cut": bool(z.imag == 0 and z.real < 0)}


# ---------------------------------------------------------------- rim collisions

def collision_points() -> list[dict]:
    """The two transverse double roots on the rim (GATE_6_PROOF.md, section 2)."""
    a = np.sqrt(2)-1
    lam = (CAP_RADIUS**-2-25/4)/a
    return [{"t": 103/115, "lam": lam, "s": 2-np.sqrt(2), "K": 1/CAP_RADIUS},
            {"t": 97/85, "lam": lam, "s": 2-np.sqrt(2), "K": -1/CAP_RADIUS}]


def rim_slice(t: float, lam: float, s: float) -> dict:
    """Roots and meridian periods of the scaled regulated denominator.

    N(K) = H*(-1+(1-t)*K) + (i/lam)*(K-5/2)**2 with H = 1-(1+Q)*s. Periods
    are coefficients of dK/(2*pi*i) for 4*Q*J**2/(K*(1-K)*N**3) at zero
    regulator in the scaled chart.
    """
    if lam <= 0:
        raise ValueError("lambda must be positive")
    H = 1-(1+ROTATION)*s
    b = 1-t
    linear, constant = 5+1j*lam*H*b, 25/4+1j*lam*H
    discriminant = linear**2-4*constant
    split = np.sqrt(discriminant)
    roots = np.array([(linear+split)/2, (linear-split)/2])

    def N(K):
        return H*(-1+b*K)+1j*(K-2.5)**2/lam

    def density(K):
        return 4*ROTATION*(-1+b*K)**2/(K*(1-K)*N(K)**3)

    paired = 4*ROTATION*(t*t/N(1)**3-1/N(0)**3)
    separation = abs(roots[0]-roots[1])
    singles = None
    clearance = min(separation, *abs(roots), *abs(roots-1))/3
    if clearance > 1e-3:
        unit = np.exp(2j*pi*np.arange(1024)/1024)
        singles = [complex(np.mean(clearance*unit*density(root+clearance*unit)))
                   for root in roots]
    scaled = abs(roots)*CAP_RADIUS
    inside = int(np.count_nonzero(scaled < 1))
    on_rim = bool(separation < .6 and abs(scaled[0]-1) < .03)
    regime = ("collision" if on_rim else ("outside", "one-meridian", "paired")[inside])
    return {
        "t": t, "lam": lam, "s": s,
        "H": [H.real, H.imag],
        "roots": [[r.real, r.imag] for r in roots],
        "scaled_modulus": scaled.tolist(),
        "separation": float(separation),
        "discriminant": float(abs(discriminant)),
        "inside": inside,
        "regime": regime,
        "paired_period": [paired.real, paired.imag],
        "single_periods": None if singles is None else [[v.real, v.imag] for v in singles],
        "radius": CAP_RADIUS,
    }


# ---------------------------------------------------------------- command line

def main(argv=None) -> None:
    parser = argparse.ArgumentParser(prog="geocusp-solver", description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="command")
    sub.add_parser("gui", help="open the desktop window (default)")

    fun = sub.add_parser("functional", help="evaluate the three functionals on a weight")
    fun.add_argument("--logs", type=float, nargs="*", default=[pi**2, 0, 1],
                     help="coefficients of log(z)**k, k=0,1,... (default: pi^2 0 1)")
    fun.add_argument("--amplitude", type=float, default=0.)
    fun.add_argument("--exponent", type=float, default=0.)
    fun.add_argument("--loop-order", type=int, default=2)
    fun.add_argument("--radius", type=float, default=1.3)
    fun.add_argument("--eta", type=float, default=.01)
    fun.add_argument("--no-edge", action="store_true",
                     help="skip the Lorentzian edge chain (it takes about a second)")

    edg = sub.add_parser("edge", help="Lorentzian edge chain at a point of a null line")
    edg.add_argument("--logs", type=float, nargs="*", default=[pi**2, 0, 1])
    edg.add_argument("--amplitude", type=float, default=0.)
    edg.add_argument("--exponent", type=float, default=0.)
    edg.add_argument("--two-loop", action="store_true", help="use the two-loop weight F2")
    edg.add_argument("--cusp", type=int, default=2)
    edg.add_argument("--edge", choices=("next", "prev"), default="next")
    edg.add_argument("--position", type=float, default=.4,
                     help="along the null line; the edge segment is 0 < position < 1")
    edg.add_argument("--region", choices=tuple(REGIONS), default="--",
                     help="signs of the two diagonals s and t")

    two = sub.add_parser("two-loop", help="the two-loop weight F2 and the three-loop coefficient")
    two.add_argument("--no-edge", action="store_true",
                     help="skip the edge chain on the continued sheets (about a minute)")

    con = sub.add_parser("contour", help="survey the explicit integration contour")
    con.add_argument("--points", type=int, default=5000)
    con.add_argument("--seed", type=int, default=0)

    cd = sub.add_parser("cd", help="integrate the CD weight at a complex z")
    cd.add_argument("z", type=complex, help="for example 0.35+0.8j")

    rim = sub.add_parser("rim", help="roots and meridian periods on one slice")
    rim.add_argument("--t", type=float, default=1.)
    rim.add_argument("--lam", type=float, default=460.)
    rim.add_argument("--s", type=float, default=.27)

    args = parser.parse_args(argv)
    if args.command in (None, "gui"):
        try:
            from .solver_gui import run
        except ImportError as error:
            parser.error(f"the desktop window needs Tk and Matplotlib ({error}); install your "
                         "platform's Tk package for Python, or use the functional, edge, "
                         "two-loop, contour, cd and rim subcommands")
        run()
        return
    try:
        if args.command == "functional":
            weight = Weight(tuple(args.logs) or (0.,), args.amplitude, args.exponent)
            result = solve_functional(weight, args.loop_order, args.radius, args.eta,
                                      edge=not args.no_edge)
        elif args.command == "edge":
            weight = ("two-loop" if args.two_loop else
                      Weight(tuple(args.logs) or (0.,), args.amplitude, args.exponent))
            result = edge_chain(weight, args.cusp, args.edge, args.position, args.region)
        elif args.command == "two-loop":
            result = solve_two_loop(edge=not args.no_edge)
        elif args.command == "contour":
            result = contour_survey(args.points, args.seed)
        elif args.command == "cd":
            result = solve_cd(args.z)
        else:
            result = rim_slice(args.t, args.lam, args.s)
    except ValueError as error:
        parser.error(str(error))
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
