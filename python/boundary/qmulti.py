"""Exact arithmetic in Q(sqrt 2, sqrt 5), the field of degree 4 that the Platonic pieces on the regular 4-polytopes need.

An element is (n0 + n1 sqrt 2 + n2 sqrt 5 + n3 sqrt 10) / d with integers n0..n3 and d > 0 in lowest terms (the gcd of all five is 1),
so equal numbers have equal fields and hash alike. Elements of the subfields Q, Q(sqrt 2) and Q(sqrt 5) are the ones with zero
coefficients: there is one class for every piece and container, and the zero parts cost only integer additions.

Signs are decided exactly. Write x = p + q sqrt 5 with p = (n0 + n1 sqrt 2) / d and q = (n2 + n3 sqrt 2) / d in Q(sqrt 2). If p and q have
opposite signs, sign(x) = sign(p) when p^2 > 5 q^2 and sign(q) otherwise; p^2 - 5 q^2 is an element of Q(sqrt 2) whose sign is decided the same
way from r^2 - 2 t^2 for r + t sqrt 2. None of these differences is ever 0 for a nonzero argument, because sqrt 2 and sqrt 5 are irrational and
independent. A binary64 evaluation decides first when it is far from 0 compared with its proven error bound (see sign()).

Mixed arithmetic. The binary operators accept int, Fraction and general.Q5 operands (a Q5 is a + b sqrt 5), so the golden containers' Q5
data and the pieces' coordinates combine without conversion code at the call sites.
"""

import math
from fractions import Fraction

_R2, _R5, _R10 = math.sqrt(2), math.sqrt(5), math.sqrt(10)  # correctly rounded binary64 square roots

_S = 10 ** 45
SQRT2_LO = Fraction(math.isqrt(2 * _S * _S), _S)
SQRT2_HI = SQRT2_LO + Fraction(1, _S)
SQRT5_LO = Fraction(math.isqrt(5 * _S * _S), _S)
SQRT5_HI = SQRT5_LO + Fraction(1, _S)
SQRT10_LO = Fraction(math.isqrt(10 * _S * _S), _S)
SQRT10_HI = SQRT10_LO + Fraction(1, _S)


def _sgn(x):
    return (x > 0) - (x < 0)


def _sign2(r, t):
    """Sign of r + t sqrt 2 for integers r, t."""
    sr, st = _sgn(r), _sgn(t)
    if st == 0:
        return sr
    if sr == 0 or sr == st:
        return st
    # opposite signs: compare r^2 with 2 t^2 (never equal for nonzero integers)
    return sr if r * r > 2 * t * t else st


class QM:
    """(n0 + n1 sqrt 2 + n2 sqrt 5 + n3 sqrt 10) / d."""

    __slots__ = ("n0", "n1", "n2", "n3", "d")

    def __init__(self, a=0, b=0, c=0, d=0):
        """a + b sqrt 2 + c sqrt 5 + d sqrt 10 with rational a, b, c, d (ints or Fractions)."""
        a, b, c, d = (x if isinstance(x, Fraction) else Fraction(x) for x in (a, b, c, d))
        den = math.lcm(a.denominator, b.denominator, c.denominator, d.denominator)
        n = [a.numerator * (den // a.denominator), b.numerator * (den // b.denominator), c.numerator * (den // c.denominator), d.numerator * (den // d.denominator)]
        g = math.gcd(n[0], n[1], n[2], n[3], den)
        self.n0, self.n1, self.n2, self.n3, self.d = n[0] // g, n[1] // g, n[2] // g, n[3] // g, den // g

    @staticmethod
    def _raw(n0, n1, n2, n3, d):
        g = math.gcd(n0, n1, n2, n3, d)
        o = object.__new__(QM)
        if g != 1:
            n0, n1, n2, n3, d = n0 // g, n1 // g, n2 // g, n3 // g, d // g
        o.n0, o.n1, o.n2, o.n3, o.d = n0, n1, n2, n3, d
        return o

    @staticmethod
    def of(x):
        if isinstance(x, QM):
            return x
        if isinstance(x, int):
            return QM._raw(x, 0, 0, 0, 1)
        if isinstance(x, Fraction):
            return QM._raw(x.numerator, 0, 0, 0, x.denominator)
        if hasattr(x, "n1") and hasattr(x, "n2") and hasattr(x, "d"):  # general.Q5: (n1 + n2 sqrt 5) / d
            return QM._raw(x.n1, 0, x.n2, 0, x.d)
        raise TypeError(f"not a number of Q(sqrt 2, sqrt 5): {x!r}")

    def coeffs(self):
        """The four rational coefficients of 1, sqrt 2, sqrt 5, sqrt 10."""
        d = self.d
        return (Fraction(self.n0, d), Fraction(self.n1, d), Fraction(self.n2, d), Fraction(self.n3, d))

    # ------------------------------------------------------------ arithmetic
    def __add__(self, o):
        if not isinstance(o, QM):
            try:
                o = QM.of(o)
            except TypeError:
                return NotImplemented
        if self.d == o.d:
            return QM._raw(self.n0 + o.n0, self.n1 + o.n1, self.n2 + o.n2, self.n3 + o.n3, self.d)
        a, b = self.d, o.d
        return QM._raw(self.n0 * b + o.n0 * a, self.n1 * b + o.n1 * a, self.n2 * b + o.n2 * a, self.n3 * b + o.n3 * a, a * b)

    __radd__ = __add__

    def __neg__(self):
        o = object.__new__(QM)
        o.n0, o.n1, o.n2, o.n3, o.d = -self.n0, -self.n1, -self.n2, -self.n3, self.d
        return o

    def __sub__(self, o):
        if not isinstance(o, QM):
            try:
                o = QM.of(o)
            except TypeError:
                return NotImplemented
        if self.d == o.d:
            return QM._raw(self.n0 - o.n0, self.n1 - o.n1, self.n2 - o.n2, self.n3 - o.n3, self.d)
        a, b = self.d, o.d
        return QM._raw(self.n0 * b - o.n0 * a, self.n1 * b - o.n1 * a, self.n2 * b - o.n2 * a, self.n3 * b - o.n3 * a, a * b)

    def __rsub__(self, o):
        try:
            return QM.of(o) - self
        except TypeError:
            return NotImplemented

    def __mul__(self, o):
        if not isinstance(o, QM):
            if isinstance(o, int):
                return QM._raw(self.n0 * o, self.n1 * o, self.n2 * o, self.n3 * o, self.d)
            try:
                o = QM.of(o)
            except TypeError:
                return NotImplemented
        a0, a1, a2, a3 = self.n0, self.n1, self.n2, self.n3
        b0, b1, b2, b3 = o.n0, o.n1, o.n2, o.n3
        return QM._raw(
            a0 * b0 + 2 * a1 * b1 + 5 * a2 * b2 + 10 * a3 * b3,
            a0 * b1 + a1 * b0 + 5 * (a2 * b3 + a3 * b2),
            a0 * b2 + a2 * b0 + 2 * (a1 * b3 + a3 * b1),
            a0 * b3 + a3 * b0 + a1 * b2 + a2 * b1,
            self.d * o.d,
        )

    __rmul__ = __mul__

    def inv(self):
        # x = p + q sqrt 5, p = a0 + a1 sqrt 2, q = a2 + a3 sqrt 2 (the common denominator d is carried along)
        a0, a1, a2, a3 = self.n0, self.n1, self.n2, self.n3
        # N = p^2 - 5 q^2 = (a0^2 + 2 a1^2 - 5 a2^2 - 10 a3^2) + (2 a0 a1 - 10 a2 a3) sqrt 2 = r + t sqrt 2
        r = a0 * a0 + 2 * a1 * a1 - 5 * a2 * a2 - 10 * a3 * a3
        t = 2 * a0 * a1 - 10 * a2 * a3
        m = r * r - 2 * t * t
        if m == 0:
            raise ZeroDivisionError("QM inverse of 0")
        # 1/x = (p - q sqrt 5) / N = (p - q sqrt 5)(r - t sqrt 2) / m
        pr = (a0 * r - 2 * a1 * t, a1 * r - a0 * t)  # p (r - t sqrt 2) = pr0 + pr1 sqrt 2
        qr = (a2 * r - 2 * a3 * t, a3 * r - a2 * t)
        # (pr0 + pr1 s2) - (qr0 + qr1 s2) s5  ->  n0 = pr0, n1 = pr1, n2 = -qr0, n3 = -qr1; the factor d: 1/x = d / (...)
        if m < 0:
            return QM._raw(-pr[0] * self.d, -pr[1] * self.d, qr[0] * self.d, qr[1] * self.d, -m)
        return QM._raw(pr[0] * self.d, pr[1] * self.d, -qr[0] * self.d, -qr[1] * self.d, m)

    def __truediv__(self, o):
        if isinstance(o, QM):
            return self * o.inv()
        try:
            return self * QM.of(o).inv()
        except TypeError:
            return NotImplemented

    def __rtruediv__(self, o):
        try:
            return QM.of(o) * self.inv()
        except TypeError:
            return NotImplemented

    # ------------------------------------------------------------ sign and order
    def sign(self):
        n0, n1, n2, n3 = self.n0, self.n1, self.n2, self.n3
        # Fast path, rigorous: each term is computed with relative error under 4e-16 (a correctly rounded conversion of an integer,
        # a correctly rounded constant, one rounded product); three additions add at most 3.4e-16 of the sum of magnitudes M. The
        # float sum is therefore within 2e-15 M of the true value, and its sign is taken when it exceeds 1e-12 M.
        try:
            t0, t1, t2, t3 = float(n0), float(n1) * _R2, float(n2) * _R5, float(n3) * _R10
            m = abs(t0) + abs(t1) + abs(t2) + abs(t3)
            f = t0 + t1 + t2 + t3
            if math.isfinite(m) and abs(f) > 1e-12 * m:
                return 1 if f > 0 else -1
        except OverflowError:
            pass
        sp, sq = _sign2(n0, n1), _sign2(n2, n3)
        if sq == 0:
            return sp
        if sp == 0 or sp == sq:
            return sq
        # opposite signs: sign(p^2 - 5 q^2) = sign of (n0^2 + 2 n1^2 - 5 n2^2 - 10 n3^2) + (2 n0 n1 - 10 n2 n3) sqrt 2
        k = _sign2(n0 * n0 + 2 * n1 * n1 - 5 * n2 * n2 - 10 * n3 * n3, 2 * n0 * n1 - 10 * n2 * n3)
        return sp if k > 0 else sq

    def is_zero(self):
        return self.n0 == 0 and self.n1 == 0 and self.n2 == 0 and self.n3 == 0

    def _cmp(self, o):
        return (self - QM.of(o)).sign()

    def __eq__(self, o):
        if isinstance(o, QM):
            return self.n0 == o.n0 and self.n1 == o.n1 and self.n2 == o.n2 and self.n3 == o.n3 and self.d == o.d
        try:
            o = QM.of(o)
        except TypeError:
            return NotImplemented
        return self == o

    def __ne__(self, o):
        r = self.__eq__(o)
        return r if r is NotImplemented else not r

    def __hash__(self):
        return hash((self.n0, self.n1, self.n2, self.n3, self.d))

    def __lt__(self, o):
        try:
            return self._cmp(o) < 0
        except TypeError:
            return NotImplemented

    def __le__(self, o):
        try:
            return self._cmp(o) <= 0
        except TypeError:
            return NotImplemented

    def __gt__(self, o):
        try:
            return self._cmp(o) > 0
        except TypeError:
            return NotImplemented

    def __ge__(self, o):
        try:
            return self._cmp(o) >= 0
        except TypeError:
            return NotImplemented

    # ------------------------------------------------------------ floats
    def bounds(self):
        """Rational (lo, hi) with lo <= value <= hi, about 1e-40 apart."""
        d = self.d
        lo = Fraction(self.n0, d)
        hi = lo
        for n, (a, b) in ((self.n1, (SQRT2_LO, SQRT2_HI)), (self.n2, (SQRT5_LO, SQRT5_HI)), (self.n3, (SQRT10_LO, SQRT10_HI))):
            if n >= 0:
                lo += Fraction(n, d) * a
                hi += Fraction(n, d) * b
            else:
                lo += Fraction(n, d) * b
                hi += Fraction(n, d) * a
        return lo, hi

    def __float__(self):
        d = self.d
        return (float(self.n0) + float(self.n1) * _R2 + float(self.n2) * _R5 + float(self.n3) * _R10) / d if max(abs(self.n0), abs(self.n1), abs(self.n2), abs(self.n3)) < 2 ** 1000 else float(sum(self.bounds()) / 2)

    def fl(self):
        """(value, err): a binary64 value and a rigorous bound on its distance from the exact number (2e-15 of the sum of the
        magnitudes of the four terms, plus the division's rounding); inf when the integers are too large for binary64."""
        try:
            d = float(self.d)
            t0, t1, t2, t3 = float(self.n0), float(self.n1) * _R2, float(self.n2) * _R5, float(self.n3) * _R10
        except OverflowError:
            return float(sum(self.bounds()) / 2), math.inf
        m = abs(t0) + abs(t1) + abs(t2) + abs(t3)
        if not math.isfinite(m):
            return float(sum(self.bounds()) / 2), math.inf
        return (t0 + t1 + t2 + t3) / d, 3e-15 * m / d

    def __repr__(self):
        return f"QM({self.n0} + {self.n1}*sqrt2 + {self.n2}*sqrt5 + {self.n3}*sqrt10)/{self.d}"

    def to_json(self):
        """[a, b, c, d] as rational strings: a + b sqrt 2 + c sqrt 5 + d sqrt 10."""
        return [str(x) for x in self.coeffs()]


SQRT2 = QM(0, 1)
SQRT5 = QM(0, 0, 1)
SQRT10 = QM(0, 0, 0, 1)
PHI = QM(Fraction(1, 2), 0, Fraction(1, 2))
