"""Exact expansion at x = -1 of combinations of harmonic polylogarithms with
letters 0 and -1, in rational arithmetic.

Purpose: to prove, not measure, that a weight F(x) given in such functions
is holomorphic at z = -1 on the principal sheet (docs/THREE_LOOP_WEIGHT.md).

With y = 1 + x the forms are dx/x = dy/(y-1) and dx/(1+x) = dy/y, and
H(a_1..a_n; x) is the iterated integral of the corresponding word from the
tangential base point at y = 1.  The path is cut in three: a half circle
round y = 1 (powers of i pi), the segment from 1 to 0 (shuffle-regularised
multiple zeta values), and the piece from 0 to y.  By the composition rule
for iterated integrals

    F = sum_u d_u I_{0 -> y}(u),

and F is holomorphic at y = 0 if d_u = 0 for every word u that ends in
dy/y.  The d_u are computed in Q[i pi, zeta3, zeta5].

Numbers: a Poly is {(a, b, c): Fraction} for (i pi)^a zeta3^b zeta5^c, with
zeta2 = pi^2/6 = -(i pi)^2/6.  Multiple zeta values up to weight 6 are
reduced to that ring by solving the double shuffle relations (shuffle
product, stuffle product, Hoffman's relation) exactly; the solver checks
that the system has full rank and is consistent, so no table is assumed.
"""
from __future__ import annotations

from fractions import Fraction
from functools import lru_cache
from itertools import product
import re


# ---------------------------------------------------------------- polynomials
class Poly:
    __slots__ = ("t",)

    def __init__(self, t=None):
        self.t = {k: v for k, v in (t or {}).items() if v != 0}

    @staticmethod
    def const(c):
        return Poly({(0, 0, 0): Fraction(c)})

    @staticmethod
    def gen(i):
        k = [0, 0, 0]
        k[i] = 1
        return Poly({tuple(k): Fraction(1)})

    def _c(self, o):
        return o if isinstance(o, Poly) else Poly.const(o)

    def __add__(self, o):
        o = self._c(o)
        t = dict(self.t)
        for k, v in o.t.items():
            t[k] = t.get(k, 0) + v
        return Poly(t)

    __radd__ = __add__

    def __neg__(self):
        return Poly({k: -v for k, v in self.t.items()})

    def __sub__(self, o):
        return self + (-self._c(o))

    def __rsub__(self, o):
        return self._c(o) - self

    def __mul__(self, o):
        o = self._c(o)
        t = {}
        for k1, v1 in self.t.items():
            for k2, v2 in o.t.items():
                k = (k1[0] + k2[0], k1[1] + k2[1], k1[2] + k2[2])
                t[k] = t.get(k, 0) + v1 * v2
        return Poly(t)

    __rmul__ = __mul__

    def __truediv__(self, o):
        return Poly({k: v / Fraction(o) for k, v in self.t.items()})

    def __pow__(self, n):
        out = Poly.const(1)
        for _ in range(n):
            out = out * self
        return out

    def is_zero(self):
        return not self.t

    def __repr__(self):
        if not self.t:
            return "0"
        names = ("(i pi)", "z3", "z5")
        parts = []
        for k, v in sorted(self.t.items()):
            mono = "*".join(f"{n}^{e}" for n, e in zip(names, k) if e)
            parts.append(f"{v}" + ("*" + mono if mono else ""))
        return " + ".join(parts)

    def numeric(self):
        """Value as an mpmath number at the working precision."""
        import mpmath as mp
        tot = mp.mpc(0)
        for (a, b, c), v in self.t.items():
            tot += (mp.mpf(v.numerator) / v.denominator) * (1j * mp.pi) ** a * mp.zeta(3) ** b * mp.zeta(5) ** c
        return tot

    def __complex__(self):
        return complex(self.numeric())


IPI, Z3, Z5 = Poly.gen(0), Poly.gen(1), Poly.gen(2)
ZETA2 = -(IPI * IPI) / 6
PI2 = -(IPI * IPI)


# ------------------------------------------------- words: 0 = du/u, 1 = du/(1-u)
def shuffle(u, v):
    """Shuffle product as {word: multiplicity}."""
    if not u:
        return {v: 1}
    if not v:
        return {u: 1}
    out = {}
    for w, m in shuffle(u[1:], v).items():
        out[(u[0],) + w] = out.get((u[0],) + w, 0) + m
    for w, m in shuffle(u, v[1:]).items():
        out[(v[0],) + w] = out.get((v[0],) + w, 0) + m
    return out


def word_of(comp):
    """zeta(n1, ..., nr), n1 the largest summation index -> word 0^{n1-1} 1 ..."""
    out = ()
    for n in comp:
        out += (0,) * (n - 1) + (1,)
    return out


def comp_of(word):
    out, run = [], 1
    for letter in word:
        if letter == 0:
            run += 1
        else:
            out.append(run)
            run = 1
    assert run == 1, "word must end in 1"
    return tuple(out)


def stuffle(a, b):
    if not a:
        return {b: 1}
    if not b:
        return {a: 1}
    out = {}
    for w, m in stuffle(a[1:], b).items():
        k = (a[0],) + w
        out[k] = out.get(k, 0) + m
    for w, m in stuffle(a, b[1:]).items():
        k = (b[0],) + w
        out[k] = out.get(k, 0) + m
    for w, m in stuffle(a[1:], b[1:]).items():
        k = (a[0] + b[0],) + w
        out[k] = out.get(k, 0) + m
    return out


def convergent(word):
    return len(word) >= 2 and word[0] == 0 and word[-1] == 1


def to_convergent(word):
    """Shuffle regularisation with zeta(0) = zeta(1) = 0: a word as a rational
    combination of convergent words (and the empty word)."""
    return dict(_to_conv(tuple(word)))


@lru_cache(maxsize=None)
def _to_conv(word):
    if not word:
        return (((), Fraction(1)),)
    if word[0] == 1:                      # leading 1's: divergent at u = 1
        p = 0
        while p < len(word) and word[p] == 1:
            p += 1
        rest = word[p:]
        if not rest:
            return ()
        # (1) sh (1^{p-1} rest) = p * 1^p rest + others with fewer leading 1's
        acc = {}
        for w, m in shuffle((1,), (1,) * (p - 1) + rest).items():
            if w == word:
                assert m == p
                continue
            for k, v in _to_conv(w):
                acc[k] = acc.get(k, 0) - Fraction(m, p) * v
        return tuple(acc.items())
    if word[-1] == 0:                     # trailing 0's: divergent at u = 0
        q = 0
        while q < len(word) and word[-1 - q] == 0:
            q += 1
        rest = word[:len(word) - q]
        if not rest:
            return ()
        acc = {}
        for w, m in shuffle((0,), rest + (0,) * (q - 1)).items():
            if w == word:
                assert m == q
                continue
            for k, v in _to_conv(w):
                acc[k] = acc.get(k, 0) - Fraction(m, q) * v
        return tuple(acc.items())
    return ((word, Fraction(1)),)


# ------------------------------------------------------------- the MZV table
def monomials(weight):
    out = []
    for a in range(0, weight + 1, 2):
        for b in range(0, weight // 3 + 1):
            for c in range(0, weight // 5 + 1):
                if a + 3 * b + 5 * c == weight:
                    out.append((a, b, c))
    return out


def solve_mzv(max_weight=6, log=None):
    """{convergent word: Poly} for all weights up to max_weight."""
    table = {(): Poly.const(1)}
    generators = {(0, 1): ZETA2, (0, 0, 1): Z3, (0, 0, 0, 0, 1): Z5}
    for n in range(2, max_weight + 1):
        words = [w for w in product((0, 1), repeat=n) if convergent(w)]
        index = {w: i for i, w in enumerate(words)}
        monos = monomials(n)
        mi = {m: i for i, m in enumerate(monos)}
        rows = []

        def add(lhs, rhs):
            row = [Fraction(0)] * len(words)
            for w, m in lhs.items():
                row[index[w]] += Fraction(m)
            r = [Fraction(0)] * len(monos)
            for k, v in rhs.t.items():
                r[mi[k]] += v
            rows.append(row + r)

        for w, val in generators.items():
            if len(w) == n:
                add({w: 1}, val)
        for k in range(2, n - 1):
            for u in product((0, 1), repeat=k):
                if not convergent(u):
                    continue
                for v in product((0, 1), repeat=n - k):
                    if not convergent(v):
                        continue
                    rhs = table[u] * table[v]
                    add(shuffle(u, v), rhs)
                    add({word_of(c): m for c, m in stuffle(comp_of(u), comp_of(v)).items()}, rhs)
        for u in product((0, 1), repeat=n - 1):       # Hoffman: zeta(1) * u - 1 sh u
            if not convergent(u):
                continue
            lhs = {}
            for c, m in stuffle((1,), comp_of(u)).items():
                w = word_of(c)
                lhs[w] = lhs.get(w, 0) + m
            for w, m in shuffle((1,), u).items():
                lhs[w] = lhs.get(w, 0) - m
            lhs = {w: m for w, m in lhs.items() if m}
            assert all(convergent(w) for w in lhs), "divergent terms must cancel"
            add(lhs, Poly())
        nw = len(words)
        # Gauss-Jordan over Q
        piv = 0
        for col in range(nw):
            r = next((i for i in range(piv, len(rows)) if rows[i][col] != 0), None)
            if r is None:
                raise ArithmeticError(f"weight {n}: relations do not determine {words[col]}")
            rows[piv], rows[r] = rows[r], rows[piv]
            f = rows[piv][col]
            rows[piv] = [x / f for x in rows[piv]]
            for i in range(len(rows)):
                if i != piv and rows[i][col] != 0:
                    g = rows[i][col]
                    rows[i] = [x - g * y for x, y in zip(rows[i], rows[piv])]
            piv += 1
        for row in rows[nw:]:
            if any(x != 0 for x in row):
                raise ArithmeticError(f"weight {n}: inconsistent relations")
        for w in words:
            row = rows[index[w]]
            table[w] = Poly({m: row[nw + j] for j, m in enumerate(monos)})
        if log is not None:
            log.append((n, nw, len(rows)))
    return table


_TABLE = None


def zeta_reg(word):
    """Shuffle-regularised zeta of a word in (0 = du/u, 1 = du/(1-u))."""
    global _TABLE
    if _TABLE is None:
        _TABLE = solve_mzv(6)
    out = Poly()
    for w, c in to_convergent(word).items():
        out = out + _TABLE[w] * c
    return out


# ------------------------------------------------------------------ parsing
def expand_index(index):
    out = []
    for m in index:
        out.extend([0] if m == 0 else [0] * (abs(m) - 1) + [1 if m > 0 else -1])
    return tuple(out)


class _Term:
    """Linear combination of H-words with Poly coefficients."""

    def __init__(self, d=None):
        self.d = d or {}

    def _lift(self, o):
        return o if isinstance(o, _Term) else _Term({None: o if isinstance(o, Poly) else Poly.const(o)})

    def __add__(self, o):
        o = self._lift(o)
        d = dict(self.d)
        for k, v in o.d.items():
            d[k] = d.get(k, Poly()) + v
        return _Term(d)

    __radd__ = __add__

    def __neg__(self):
        return _Term({k: -v for k, v in self.d.items()})

    def __sub__(self, o):
        return self + (-self._lift(o))

    def __rsub__(self, o):
        return self._lift(o) - self

    def __mul__(self, o):
        o = self._lift(o)
        d = {}
        for k1, v1 in self.d.items():
            for k2, v2 in o.d.items():
                if k1 is not None and k2 is not None:
                    raise ValueError("products of harmonic polylogarithms are not supported")
                k = k1 if k1 is not None else k2
                d[k] = d.get(k, Poly()) + v1 * v2
        return _Term(d)

    __rmul__ = __mul__

    def __truediv__(self, o):
        if isinstance(o, _Term):
            if set(o.d) != {None} or set(o.d[None].t) != {(0, 0, 0)}:
                raise ValueError("division by a non-rational quantity")
            o = o.d[None].t[(0, 0, 0)]
        return _Term({k: v / o for k, v in self.d.items()})

    def __pow__(self, n):
        out = _Term({None: Poly.const(1)})
        for _ in range(n):
            out = out * self
        return out


def parse_exact(text):
    """Mathematica expression in HPL[{...}, x] -> {word in (0,-1): Poly}; the
    empty word carries the constant."""
    words = []

    def replace(match):
        word = expand_index(tuple(int(t) for t in match.group(1).split(",")))
        if any(letter not in (0, -1) for letter in word):
            raise ValueError("only the letters 0 and -1 are supported")
        words.append(word)
        return "H[w%dw]" % (len(words) - 1)

    source = re.sub(r"HPL\[\{([^}]*)\},\s*x\]", replace, text.replace("\n", " "))
    source = source.replace("Zeta[3]", "Z3").replace("Zeta[5]", "Z5")
    source = re.sub(r"Pi\^(\d+)", lambda m: "PIPOW%sE" % m.group(1), source).replace("Pi", "PIPOW1E")
    if re.search(r"\d\.|\.\d|\d[eE][-+]?\d", source):
        raise ValueError("only integers and their ratios are supported")
    source = re.sub(r"(Z[35])\^(\d+)", lambda m: "(" + "*".join([m.group(1)] * int(m.group(2))) + ")", source)
    if "^" in source:
        raise ValueError("unsupported power")
    source = re.sub(r"(?<![\w.])(\d+)(?![\w.])", r"I(\1)", source)
    source = re.sub(r"PIPOW(\d+)E", r"PI(\1)", source)
    source = re.sub(r"H\[w(\d+)w\]", r"H(\1)", source)
    if re.search(r"[\[\]{}]", source):
        raise ValueError("unrecognised function in the expression")
    one = _Term({None: Poly.const(1)})
    env = {"H": lambda i: _Term({words[i]: Poly.const(1)}), "Z3": one * Z3, "Z5": one * Z5, "PI": _pi_power,
           "I": lambda n: one * Poly.const(n),
           "__builtins__": {}}
    value = eval(source, env)                       # noqa: S307 - the source is checked above
    out = {}
    for k, v in value.d.items():
        out[() if k is None else k] = v
    return out


def _pi_power(n):
    """Only even powers of pi lie in the ring."""
    if n % 2:
        raise ValueError("odd power of pi")
    return _Term({None: PI2 ** (n // 2)})


# ------------------------------------------------- expansion at x = -1 (y = 0)
def transport(v, sign=+1):
    """T(v) = sum over v = v2 v1 of I_{1 -> 0}(v2) I_{half circle}(v1) for a
    word v in the forms (0: dy/y, 1: dy/(y-1)), outermost first.  sign = +1
    for the lip Im L = +pi."""
    out = Poly()
    n = len(v)
    for k in range(n + 1):
        v2, v1 = v[:k], v[k:]
        if any(letter != 1 for letter in v1):
            continue
        m = len(v1)
        half = (IPI * sign) ** m
        for j in range(2, m + 1):
            half = half / j
        # y = 1 - u: dy/y = -du/(1-u), dy/(y-1) = du/u
        zword = tuple(1 if letter == 0 else 0 for letter in v2)
        segment = zeta_reg(zword) * (-1) ** sum(1 for letter in v2 if letter == 0)
        out = out + segment * half
    return out


def expansion_at_minus_one(coefficients, sign=+1):
    """{word u in (0: dy/y, 1: dy/(y-1)): Poly}, the coefficients d_u of
    F = sum_u d_u I_{0 -> y}(u), y = 1 + x."""
    d = {}
    cache = {}
    for xword, c in coefficients.items():
        w = tuple(1 if a == 0 else 0 for a in xword)        # dx/x = dy/(y-1), dx/(1+x) = dy/y
        for j in range(len(w) + 1):
            u, v = w[:j], w[j:]
            if v not in cache:
                cache[v] = transport(v, sign)
            d[u] = d.get(u, Poly()) + c * cache[v]
    return {u: p for u, p in d.items()}


def holomorphy_report(coefficients, sign=+1):
    """value: F at x = -1 if the list `logarithmic` is empty; logarithmic: the
    non-zero d_u with u ending in dy/y (F is holomorphic at x = -1 on the
    lip Im L = sign * pi if there are none); slope: dF/dy at y = 0 in that
    case, which is the limit of F/(1 + x) when the value is zero."""
    d = expansion_at_minus_one(coefficients, sign)
    logs = {u: p for u, p in d.items() if u and u[-1] == 0 and not p.is_zero()}
    slope = Poly()
    for u, p in d.items():                      # I(0^k 1) = -y + O(y^2); two or more 1's are O(y^2)
        if u and u[-1] == 1 and all(letter == 0 for letter in u[:-1]):
            slope = slope - p
    return {"value": d.get((), Poly()), "slope": slope, "logarithmic": logs,
            "conditions": sum(1 for u in d if u and u[-1] == 0), "coefficients": d}


def series_at_minus_one(d, y, terms=80):
    """sum_u d_u I_{0 -> y}(u) for |y| < 1 from the power series, for a set of
    coefficients with no logarithmic term.  Independent of the HPL evaluator."""
    import mpmath as mp

    @lru_cache(maxsize=None)
    def coeffs(u):
        if not u:
            return (mp.mpf(1),) + (mp.mpf(0),) * terms
        inner = coeffs(u[1:])
        out = [mp.mpf(0)] * (terms + 1)
        if u[0] == 0:                           # int dy/y
            if inner[0] != 0:
                raise ValueError("logarithmic term")
            for n in range(1, terms + 1):
                out[n] = inner[n] / n
        else:                                   # int dy/(y-1) = -int sum_k y^k
            run = mp.mpf(0)
            for n in range(1, terms + 1):
                run += inner[n - 1]
                out[n] = -run / n
        return tuple(out)

    total = mp.mpc(0)
    for u, p in d.items():
        if p.is_zero():
            continue
        c = coeffs(u)
        total += p.numeric() * mp.polyval(list(reversed(c)), y)
    return total
