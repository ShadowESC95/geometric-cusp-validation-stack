"""Harmonic polylogarithms with letters 0 and -1, and weights written in them.

Used to evaluate the functions F^(L)(x) of Henn, Korchemsky and Mistlberger
(arXiv:1911.10174), which are given as Mathematica expressions in
HPL[{...}, x]. The one- and two-loop expressions are short and are included
here; the three-loop files are read from the paper's ancillary directory,
which is not distributed with this repository:

    python -m geocusp.hpl /path/to/anc

Normalisation: F-bar^(L) = 4^(L+1) F^(L)_HKM, and
Gamma^(L) = -(4/L) M[F-bar^(L-1)] in powers of g^2.
"""
from __future__ import annotations

import re

import mpmath as mp
import numpy as np
from scipy.integrate import solve_ivp

F1_HKM = "Pi^2/16 + HPL[{0, 0}, x]/8"
F2_HKM = ("(-107*Pi^4)/5760 - (Pi^2*HPL[{-1, -1}, x])/16 + (Pi^2*HPL[{-1, 0}, x])/32 + "
          "(Pi^2*HPL[{0, -1}, x])/32 - (Pi^2*HPL[{0, 0}, x])/16 - HPL[{-1, -1, 0, 0}, x]/8 + "
          "HPL[{-1, 0, 0, 0}, x]/16 + HPL[{0, -1, 0, 0}, x]/16 - (3*HPL[{0, 0, 0, 0}, x])/16 - "
          "(HPL[{-1}, x]*Zeta[3])/8 + (HPL[{0}, x]*Zeta[3])/16")


def expand_index(index):
    """Compressed HPL notation: m with |m| > 1 stands for |m|-1 zeros and sign(m)."""
    out = []
    for m in index:
        out.extend([0] if m == 0 else [0] * (abs(m) - 1) + [1 if m > 0 else -1])
    return tuple(out)


def parse(text):
    """A Mathematica expression in HPL[{...}, x] -> (python source, set of words)."""
    words = set()

    def replace(match):
        word = expand_index(tuple(int(t) for t in match.group(1).split(",")))
        if any(letter not in (0, -1) for letter in word):
            raise ValueError("only the letters 0 and -1 are supported")
        words.add(word)
        return "H[%r]" % (word,)

    source = re.sub(r"HPL\[\{([^}]*)\},\s*x\]", replace, text.replace("\n", " "))
    source = source.replace("Zeta[3]", "Z3").replace("Zeta[5]", "Z5").replace("Pi", "PI").replace("^", "**")
    if re.search(r"[A-Za-z_]+\[", source.replace("H[", "")):
        raise ValueError("unrecognised function in the expression")
    return source, words


def suffix_closure(words):
    out = set()
    for word in words:
        for i in range(len(word)):
            out.add(word[i:])
    return sorted(out, key=lambda w: (len(w), w))


def series_values(words, x, terms=70):
    """H_w(x) for |x| < 1: power series in x with powers of log x.

    `words` must be closed under removing the first letter.
    H(a, w; x) = int_0^x f_a(t) H(w; t) dt, f_0 = 1/t, f_{-1} = 1/(1+t), and
    H(0,...,0; x) = log^n(x)/n!.
    """
    with mp.workdps(22):
        x = mp.mpmathify(x)
        log_x = mp.log(x)
        series = {(): {0: [mp.mpf(1)] + [mp.mpf(0)] * terms}}

        def integrate(source, letter):
            result = {}
            for k, coefficients in source.items():
                if letter == 0:
                    g, shift = coefficients, 0
                else:
                    g, running = [mp.mpf(0)] * (terms + 1), mp.mpf(0)
                    for n in range(terms + 1):
                        running = coefficients[n] - running
                        g[n] = running
                    shift = 1
                for n in range(terms + 1):
                    c, power = g[n], n + shift
                    if c == 0:
                        continue
                    if power == 0:
                        result.setdefault(k + 1, [mp.mpf(0)] * (terms + 1))[0] += c / (k + 1)
                    elif power <= terms:
                        falling = mp.mpf(1)
                        for j in range(k + 1):
                            result.setdefault(k - j, [mp.mpf(0)] * (terms + 1))[power] += (
                                c * (-1) ** j * falling / mp.mpf(power) ** (j + 1))
                            falling *= k - j
            return result

        values = {}
        for word in words:
            if word not in series:
                series[word] = integrate(series[word[1:]], word[0])
            total = mp.mpf(0)
            for k, coefficients in series[word].items():
                total += log_x ** k * mp.polyval(coefficients[::-1], x)
            values[word] = complex(total)
        return values


class HPLWeights:
    """A set of named expressions in HPLs of x, evaluated together."""

    def __init__(self, expressions: dict):
        self.code, words = {}, set()
        for name, text in expressions.items():
            source, used = parse(text)
            self.code[name] = compile(source, name, "eval")
            words |= used
        self.words = suffix_closure(words)
        self.index = {w: i for i, w in enumerate(self.words)}
        self._tail = np.array([self.index[w[1:]] if len(w) > 1 else -1 for w in self.words])
        self._letter = np.array([w[0] for w in self.words])
        self._constants = {"PI": np.pi, "Z3": float(mp.zeta(3)), "Z5": float(mp.zeta(5))}

    def _rhs(self, x, y):
        previous = np.where(self._tail >= 0, y[self._tail], 1.0)
        return np.where(self._letter == 0, 1 / x, 1 / (1 + x)) * previous

    def transport(self, y, path, derivative, rtol=1e-12):
        """Continue all HPLs along x = path(t), 0 <= t <= 1."""
        solution = solve_ivp(lambda t, v: self._rhs(path(t), v) * derivative(t), (0.0, 1.0), y,
                             method="DOP853", rtol=rtol, atol=1e-14)
        if not solution.success:
            raise ArithmeticError("continuation of the HPLs failed: " + solution.message)
        return solution.y[:, -1]

    def at(self, x):
        """All HPLs at x on the principal sheet (x not on the negative axis).

        Series at 0.4 x/|x|, then a path that is straight in log(1+x), on
        which the letter -1 is regular; it stays in the same half plane.
        """
        x = complex(x)
        start = 0.4 * x / abs(x)
        values = series_values(self.words, start)
        y = np.array([values[w] for w in self.words], dtype=complex)
        rate = np.log((1 + x) / (1 + start))
        return self.transport(y, lambda t: (1 + start) * np.exp(rate * t) - 1,
                              lambda t: (1 + start) * np.exp(rate * t) * rate)

    def around(self, y, x, centre, turns=1):
        """Continue from x once (or `turns` times) around `centre`, back to x."""
        radius = x - centre
        return self.transport(y, lambda t: centre + radius * np.exp(2j * np.pi * turns * t),
                              lambda t: radius * 2j * np.pi * turns * np.exp(2j * np.pi * turns * t))

    def value(self, name, x, y):
        H = {w: y[i] for w, i in self.index.items()}
        return complex(eval(self.code[name], {"H": H, "x": x, **self._constants}))

    def circle_mean(self, combine=None, nodes=40):
        """M[F] = (1/pi) Re int_0^pi F(e^{i phi}) d phi on the principal sheet.

        `combine(values, x)` may build further functions from the named ones.
        The substitution phi = pi (1 - u^4) removes the logarithms at phi = pi.
        """
        u, weights = np.polynomial.legendre.leggauss(nodes)
        u, weights = (u + 1) / 2, weights / 2
        totals: dict = {}
        for ui, wi in zip(u, weights):
            x = np.exp(1j * np.pi * (1 - ui ** 4))
            y = self.at(x)
            values = {name: self.value(name, x, y) for name in self.code}
            if combine is not None:
                values.update(combine(values, x))
            for name, v in values.items():
                totals[name] = totals.get(name, 0.0) + wi * 4 * ui ** 3 * v.real
        return totals


def series_values_mp(words, x, terms=130):
    """`series_values` at the current mpmath precision, returning mpmath numbers."""
    x = mp.mpmathify(x)
    log_x = mp.log(x)
    series, values = {(): {0: [mp.mpf(1)] + [mp.mpf(0)] * terms}}, {}
    for word in words:
        result = {}
        for k, coefficients in series[word[1:]].items():
            if word[0] == 0:
                g, shift = coefficients, 0
            else:
                g, running = [mp.mpf(0)] * (terms + 1), mp.mpf(0)
                for n in range(terms + 1):
                    running = coefficients[n] - running
                    g[n] = running
                shift = 1
            for n in range(terms + 1):
                c, power = g[n], n + shift
                if c == 0:
                    continue
                if power == 0:
                    result.setdefault(k + 1, [mp.mpf(0)] * (terms + 1))[0] += c / (k + 1)
                elif power <= terms:
                    falling = mp.mpf(1)
                    for j in range(k + 1):
                        result.setdefault(k - j, [mp.mpf(0)] * (terms + 1))[power] += (
                            c * (-1) ** j * falling / mp.mpf(power) ** (j + 1))
                        falling *= k - j
        series[word] = result
        values[word] = sum(log_x ** k * mp.polyval(co[::-1], x) for k, co in result.items())
    return values


def taylor_step(words, values, a, b, terms=75):
    """Continue all HPLs from a to b by their Taylor series at a.

    Needs |b - a| well inside the distance from a to 0 and to -1.
    """
    g0 = [(-1) ** m / a ** (m + 1) for m in range(terms)]
    g1 = [(-1) ** m / (1 + a) ** (m + 1) for m in range(terms)]
    coefficients, new = {(): [mp.mpf(1)] + [mp.mpf(0)] * terms}, {}
    for word in words:
        previous, g = coefficients[word[1:]], (g0 if word[0] == 0 else g1)
        c = [values[word]] + [mp.mpf(0)] * terms
        for n in range(terms):
            c[n + 1] = mp.fsum(g[m] * previous[n - m] for m in range(n + 1)) / (n + 1)
        coefficients[word] = c
        new[word] = mp.polyval(c[::-1], b - a)
    return new


def monodromy_at_minus_one(expressions: dict, combine=None, dps=45, steps=40):
    """|F after one turn around x = -1 minus F|, at high precision.

    Starts from the series at x = 0.4 i, walks in the upper half plane to
    x = -1 + 0.3 e^{0.4 i} on the principal sheet, and goes once around -1.
    Returns the differences, which vanish if F has no branch point there.
    """
    with mp.workdps(dps):
        code, words = {}, set()
        for name, text in expressions.items():
            source, used = parse(text)
            code[name] = compile(source, name, "eval")
            words |= used
        words = suffix_closure(words)

        def functions(x, values):
            env = {"H": values, "PI": mp.pi, "Z3": mp.zeta(3), "Z5": mp.zeta(5), "x": x}
            out = {name: eval(c, env) for name, c in code.items()}
            if combine is not None:
                out.update(combine(out, x))
            return out

        a = mp.mpc(0, "0.4")
        values = series_values_mp(words, a)
        target = -1 + mp.mpf("0.3") * mp.expj(mp.mpf("0.4"))
        while True:
            reach = mp.mpf("0.22") * min(abs(a), abs(1 + a))
            if abs(target - a) <= reach:
                values = taylor_step(words, values, a, target)
                break
            nxt = a + (target - a) / abs(target - a) * reach
            values, a = taylor_step(words, values, a, nxt), nxt
        before = functions(target, values)
        position = target
        for j in range(1, steps + 1):
            nxt = -1 + (target + 1) * mp.expj(2 * mp.pi * j / steps)
            values, position = taylor_step(words, values, position, nxt), nxt
        after = functions(target, values)
        return {name: float(abs(after[name] - before[name])) for name in before}


def three_loop_report(directory, nodes=40):
    """Circle means of the HKM functions against the cusp anomalous dimension,
    and the hypotheses of the double-pole theorem for the three-loop weight."""
    import os

    read = lambda name: open(os.path.join(directory, name)).read()
    weights = HPLWeights({"F1": F1_HKM, "F2": F2_HKM, "F3_planar": read("HKM_F3_P_HPL.txt"),
                          "NP1": read("HKM_F3_NP1_HPL.txt"), "NP2": read("HKM_F3_NP2_HPL.txt")})
    nonplanar = lambda v, x: {"F3_nonplanar": v["NP1"] / (1 + x) + v["NP2"]}
    means = weights.circle_mean(nonplanar, nodes)
    z2, z3, z4, z6 = (float(mp.zeta(k)) for k in (2, 3, 4, 6))
    rows = [("two loops", -2 * means["F1"], -z2 / 2),
            ("three loops", -4 / 3 * means["F2"], 11 * z4 / 8),
            ("four loops, planar", -means["F3_planar"], -(z3 ** 2 / 8 + 219 * z6 / 64)),
            ("four loops, non-planar", -means["F3_nonplanar"], -1.5 * (1.5 * z3 ** 2 + 31 * z6 / 8))]
    x = -1 + 0.3 * np.exp(0.4j)
    y = weights.at(x)

    def all_values(v):
        out = {name: weights.value(name, x, v) for name in ("F3_planar", "NP1", "NP2")}
        out["F3_nonplanar"] = out["NP1"] / (1 + x) + out["NP2"]
        return out

    base, once = all_values(y), all_values(weights.around(y, x, -1))
    monodromy = {k: abs(once[k] - base[k]) for k in ("F3_planar", "F3_nonplanar")}
    near = -1 + 1e-6 * np.exp(0.9j)
    yn = weights.at(near)
    limit = {"F3_planar": weights.value("F3_planar", near, yn),
             "NP1": weights.value("NP1", near, yn),
             "F3_nonplanar": weights.value("NP1", near, yn) / (1 + near) + weights.value("NP2", near, yn)}
    return {"cusp": rows, "monodromy_at_minus_one": monodromy, "values_near_minus_one": limit}


if __name__ == "__main__":
    import sys

    report = three_loop_report(sys.argv[1])
    for name, got, expected in report["cusp"]:
        print(f"{name:24s} from the circle mean {got:+.13f}   known {expected:+.13f}   difference {got - expected:+.1e}")
    print("monodromy around x = -1 on the principal sheet:", report["monodromy_at_minus_one"])
    print("values near x = -1:", report["values_near_minus_one"])
