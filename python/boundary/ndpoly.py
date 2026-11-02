"""Exact certifier for balls on the boundary of the d-simplex, d-orthoplex and d-demicube, d = 5..10 ("4D on 5D" to
"9D on 10D"), with developments along chains of cells in R^(d-1). BOUNDARY.md section 4D.

Problem boundary/<d><family>/ball/<P> (family simplex, orthoplex or demicube): P balls of diameter 1 on the boundary of the
polytope of edge s; the score is s. The conventions are those of packages/boundary/src/nd/container.ts and
docs/spaces/BOUNDARY.md 4D.1: cell ids, a cell's centre and its frame (Gram-Schmidt, run twice, on the projections of the standard
basis onto the cell's tangent space), a configuration {s, pieces: [{cell, u}]} with u the centre in the cell's local frame at s.

Exactness. Everything is decided in exact arithmetic: rationals for the simplex and orthoplex, and Q(sqrt k) for the
demicube (k the squarefree part of d; the simplex facets and the demicube facets meet at an angle whose cosine is 1/sqrt d). The polytope is the
BASE polytope of edge sqrt 2 (simplex: the vertices e_i of R^(d+1); orthoplex: +-e_i; demicube: (+-1/2)^d with an even number of
minus signs), with outward normals N_c of one common squared length kappa (the facet normals 1 - (d+1) e_i, sigma, sqrt d e_i and w), so no square
root of a norm is ever taken. A ball of diameter 1 on the polytope of edge s has diameter e0 / s = sqrt 2 / s in base units, so every
test compares a squared base length with 2 / s^2 times a rational: no square root appears. Only the
reported margins are converted to numbers.

Developments (4D.2). The cell c' glued to c across a ridge is unfolded into the hyperplane of c by the reflection
    sigma(x) = x - 2 (m . (x - q)) / (m . m) m,    m = N_c + N_c',   q a point of the ridge,
which fixes the ridge, maps the hyperplane of c' onto that of c and lays c' on the far side of the ridge. A node of a chain
stores a = the source as seen in its own cell's hyperplane (a child's a is sigma(a)), and every window is carried back along
the chain by the linear parts of the reflections. For a ball centre u of a cell, the chains are kept by the same sound
rule as the engine's: the ridge of neighbour g is crossed only if S_g = h_g - N_g . a > 0 and the distance from a to the ridge, at least
sqrt(S_g^2 / nu_g^2 + max_k delta_k^2) (the distance to the ridge's flat inside the hyperplane and, inside that flat, to the half-space of the cell's
other facet k), is at most the reach (half a ball's diameter, plus the strict-mode slack). A path of length L has a midpoint at
distance L / 2 from each end, so pairs and loops are found by JOINS of two half-chains ending in one cell; the hinge
(the (d-3)-face) distance by an exact nearest-point search inside the window.

The centre of a ball is read exactly from the engine's numbers: the ambient point s * centre + sum u_k frame_k (floats, as the engine computes
it) is written as an affine combination of n + 1 independent vertices of the cell with dyadic weights (rounded to 2^-60, summing to 1), and
the exact point is the same combination of the exact vertices. That point lies in the cell's hyperplane exactly and differs from the engine's by
about 1e-16 s; a point that falls outside its closed cell by less than 1e-12 s is moved onto the facet planes it violates before the exact test.
Strict validity asks for room for that difference: a pair distance or a loop must exceed 1 by more than READING_PAIR s, a clearance 1/2 by more
than READING_CLEAR s.

Usage: imported by certify.py. python3 python/boundary/ndpoly.py boundary/<d><family>/ball/<P> <config.json> certifies one configuration.
"""

import itertools
import math
import re
import sys
from fractions import Fraction

NAME_RE = re.compile(r"^(5|6|7|8|9|10)(simplex|orthoplex|demicube)$")
ID_RE = re.compile(r"^boundary/((?:5|6|7|8|9|10)(?:simplex|orthoplex|demicube))/(ball)/([1-9][0-9]{0,3})\Z")
NAMES = tuple(f"{d}{f}" for f in ("simplex", "orthoplex", "demicube") for d in (5, 6, 7, 8, 9, 10))
MAX_P = 24
MAX_S = Fraction(10) ** 4
VERSION = "boundary-nd 1"
CHECKER = "python/boundary/ndpoly.py"
# pi to 50 decimals, truncated
PI_LO = Fraction("3.14159265358979323846264338327950288419716939937510")
# reading slack, relative to s (BOUNDARY.md 4D.4); set from the measured discrepancy
READING_PAIR = Fraction(2, 10 ** 13)
READING_CLEAR = Fraction(1, 10 ** 13)
SNAP = Fraction(1, 10 ** 12)  # a centre outside its cell by less than this fraction of s is moved onto the facet planes it violates
MAX_NODES = 20_000
MAX_STEPS = 2_000_000


class NDError(Exception):
    """Malformed input or a configuration outside this certifier's scope or work limit."""


class WorkLimit(NDError):
    pass


# ------------------------------------------------------------------ numbers: Fractions, and Q(sqrt k) for the demicube


def sqrt_floor(n):
    return math.isqrt(n)


def squarefree(d):
    """d = c^2 k with k squarefree: returns (c, k)."""
    c = 1
    k = d
    f = 2
    while f * f <= k:
        while k % (f * f) == 0:
            k //= f * f
            c *= f
        f += 1
    return c, k


def make_qf(k):
    """The class of numbers (n + m sqrt k) / q with integers n, m and q > 0 (k > 1 squarefree): the value a + b sqrt k with rational a = n/q, b = m/q.
    Integer arithmetic with the denominator carried along (and reduced when the numbers grow) is several times faster than pairs of Fractions."""

    class QF:
        __slots__ = ("n", "m", "q")
        K = k

        def __init__(self, a=0, b=0, _raw=None):
            if _raw is not None:
                self.n, self.m, self.q = _raw
                return
            a = a if isinstance(a, Fraction) else Fraction(a)
            b = b if isinstance(b, Fraction) else Fraction(b)
            q = a.denominator * b.denominator // math.gcd(a.denominator, b.denominator)
            self.n, self.m, self.q = a.numerator * (q // a.denominator), b.numerator * (q // b.denominator), q

        @property
        def a(self):
            return Fraction(self.n, self.q)

        @property
        def b(self):
            return Fraction(self.m, self.q)

        @staticmethod
        def of(x):
            if isinstance(x, QF):
                return x
            if isinstance(x, int):
                return QF(_raw=(x, 0, 1))
            x = x if isinstance(x, Fraction) else Fraction(x)
            return QF(_raw=(x.numerator, 0, x.denominator))

        @staticmethod
        def _mk(n, m, q):
            if q.bit_length() > 96 or abs(n).bit_length() > 96 or abs(m).bit_length() > 96:
                g = math.gcd(n, m, q)
                if g > 1:
                    n, m, q = n // g, m // g, q // g
            return QF(_raw=(n, m, q))

        def __add__(self, o):
            if not isinstance(o, QF):
                o = QF.of(o)
            if self.q == o.q:
                return QF._mk(self.n + o.n, self.m + o.m, self.q)
            return QF._mk(self.n * o.q + o.n * self.q, self.m * o.q + o.m * self.q, self.q * o.q)

        __radd__ = __add__

        def __neg__(self):
            return QF(_raw=(-self.n, -self.m, self.q))

        def __sub__(self, o):
            if not isinstance(o, QF):
                o = QF.of(o)
            if self.q == o.q:
                return QF._mk(self.n - o.n, self.m - o.m, self.q)
            return QF._mk(self.n * o.q - o.n * self.q, self.m * o.q - o.m * self.q, self.q * o.q)

        def __rsub__(self, o):
            return QF.of(o) - self

        def __mul__(self, o):
            if not isinstance(o, QF):
                o = QF.of(o)
            return QF._mk(self.n * o.n + self.m * o.m * k, self.n * o.m + self.m * o.n, self.q * o.q)

        __rmul__ = __mul__

        def inv(self):
            d = self.n * self.n - self.m * self.m * k
            if d == 0:
                raise ZeroDivisionError("division by zero in Q(sqrt k)")
            n, m, q = self.q * self.n, -self.q * self.m, d
            if q < 0:
                n, m, q = -n, -m, -q
            return QF._mk(n, m, q)

        def __truediv__(self, o):
            return self * QF.of(o).inv()

        def __rtruediv__(self, o):
            return QF.of(o) * self.inv()

        def sign(self):
            a, b = self.n, self.m
            if b == 0:
                return (a > 0) - (a < 0)
            if a == 0:
                return (b > 0) - (b < 0)
            if (a > 0) == (b > 0):
                return 1 if a > 0 else -1
            t = a * a - b * b * k
            s = (t > 0) - (t < 0)
            return s if a > 0 else -s

        def __eq__(self, o):
            if not isinstance(o, QF):
                o = QF.of(o)
            return self.n * o.q == o.n * self.q and self.m * o.q == o.m * self.q

        def __hash__(self):
            g = math.gcd(self.n, self.m, self.q)
            return hash((self.n // g, self.m // g, self.q // g))

        def __lt__(self, o):
            return (self - o).sign() < 0

        def __le__(self, o):
            return (self - o).sign() <= 0

        def __gt__(self, o):
            return (self - o).sign() > 0

        def __ge__(self, o):
            return (self - o).sign() >= 0

        def __bool__(self):
            return self.n != 0 or self.m != 0

        def __float__(self):
            lo, hi = self.bounds(30)
            return float((lo + hi) / 2)

        def bounds(self, digits=60):
            """Rational lower and upper bounds."""
            r = Fraction(sqrt_floor(k * 10 ** (2 * digits)), 10 ** digits)
            r_hi = r + Fraction(1, 10 ** digits)
            a, b = self.a, self.b
            lo = a + b * (r if b > 0 else r_hi)
            hi = a + b * (r_hi if b > 0 else r)
            return lo, hi

    return QF


def sg(x):
    """-1, 0 or 1: the exact sign of a Fraction, an int or a Q(sqrt k) number."""
    if hasattr(x, "sign"):
        return x.sign()
    return (x > 0) - (x < 0)


def dot(u, v):
    s = 0
    for a, b in zip(u, v):
        if a and b:  # truthiness: no equality test with a number (a Q(sqrt k) test would build a number)
            s = s + a * b
    return s


def to_float(x):
    if hasattr(x, "bounds"):
        lo, hi = x.bounds(30)
        return float((lo + hi) / 2)
    return float(x)


# ------------------------------------------------------------------ the containers


def popcount(x):
    return bin(x).count("1")


class Cell:
    __slots__ = ("c", "kind", "N", "h", "nbr", "Nn", "hn", "rp", "nu2", "hinges", "Nc_dot")


class Cont:
    """The base polytope of edge sqrt 2 with unnormalised outward normals of common squared length kappa; cells built on demand."""

    def __init__(self, name):
        m = NAME_RE.match(name)
        if not m:
            raise NDError(f"unknown container {name!r}")
        self.name = name
        self.d = int(m.group(1))
        self.family = m.group(2)
        d = self.d
        self.n = d - 1
        self.D = d + 1 if self.family == "simplex" else d
        self.C = d + 1 if self.family == "simplex" else (2 ** d if self.family == "orthoplex" else 2 * d + 2 ** (d - 1))
        self.kappa = {"simplex": d * (d + 1), "orthoplex": d, "demicube": d}[self.family]
        self.sqd = None
        if self.family == "demicube":
            c, k = squarefree(d)
            if k == 1:
                self.sqd = Fraction(c)
            else:
                QF = make_qf(k)
                self.sqd = QF(0, c)
        self.odd = [m_ for m_ in range(2 ** d) if popcount(m_) % 2 == 1]
        self.odd_index = {m_: i for i, m_ in enumerate(self.odd)}
        self._cells = {}
        self._verts = {}
        self._mirrors = {}

    # ---- facets
    def signs(self, mask):
        return [-1 if (mask >> i) & 1 else 1 for i in range(self.d)]

    def normal(self, c):
        d, D = self.d, self.D
        if self.family == "simplex":
            N = [Fraction(1)] * D
            N[c] = Fraction(1 - (d + 1))
            return N, Fraction(1), 0
        if self.family == "orthoplex":
            return [Fraction(-1 if (c >> i) & 1 else 1) for i in range(d)], Fraction(1), 0
        if c < 2 * d:
            N = [Fraction(0)] * d
            N[c >> 1] = self.sqd * (1 if c & 1 else -1)
            return N, self.sqd / 2, 1
        w = self.signs(self.odd[c - 2 * d])
        return [Fraction(x) for x in w], Fraction(d - 2, 2), 2

    def neighbours(self, c):
        d = self.d
        if self.family == "simplex":
            return [g for g in range(d + 1) if g != c]
        if self.family == "orthoplex":
            return [c ^ (1 << k) for k in range(d)]
        if c < 2 * d:
            i = c >> 1
            neg = not (c & 1)
            out = []
            for j in range(d):
                if j != i:
                    out += [2 * j, 2 * j + 1]
            for k, m_ in enumerate(self.odd):
                if ((m_ >> i) & 1) == (1 if neg else 0):
                    out.append(2 * d + k)
            return out
        m_ = self.odd[c - 2 * d]
        return [2 * i + (0 if (m_ >> i) & 1 else 1) for i in range(d)]

    def ridge_point(self, c, g):
        d, D = self.d, self.D
        x = [Fraction(0)] * D
        if self.family == "simplex":
            x = [Fraction(1, d - 1)] * D
            x[c] = Fraction(0)
            x[g] = Fraction(0)
        elif self.family == "orthoplex":
            k = (c ^ g).bit_length() - 1
            for i in range(d):
                x[i] = Fraction(0) if i == k else Fraction(-1 if (c >> i) & 1 else 1, d - 1)
        elif c < 2 * d and g < 2 * d:
            x[c >> 1] = Fraction(1, 2) if c & 1 else Fraction(-1, 2)
            x[g >> 1] = Fraction(1, 2) if g & 1 else Fraction(-1, 2)
        else:
            dem = c if c < 2 * d else g
            sim = g if c < 2 * d else c
            i = dem >> 1
            w = self.signs(self.odd[sim - 2 * d])
            for k in range(d):
                x[k] = Fraction(w[k], 2) - (Fraction(0) if k == i else Fraction(w[k], d - 1))
        return x

    def cell(self, c):
        cl = self._cells.get(c)
        if cl is not None:
            return cl
        if not 0 <= c < self.C:
            raise NDError(f"cell {c} is out of range")
        cl = Cell()
        cl.c = c
        cl.N, cl.h, cl.kind = self.normal(c)
        cl.nbr = self.neighbours(c)
        F = len(cl.nbr)
        cl.Nn, cl.hn, cl.rp, cl.nu2, cl.Nc_dot = [], [], [], [], []
        for g in cl.nbr:
            Ng, hg, _ = self.normal(g)
            cl.Nn.append(Ng)
            cl.hn.append(hg)
            cl.rp.append(self.ridge_point(c, g))
            cd = dot(Ng, cl.N)
            cl.Nc_dot.append(cd)
            cl.nu2.append(self.kappa - cd * cd / self.kappa)
        d = self.d
        pairs = []
        if cl.kind == 1:
            nd = 2 * (d - 1)
            for k in range(nd):
                for l in range(k + 1, nd):
                    if d >= 6 and (cl.nbr[k] >> 1) != (cl.nbr[l] >> 1):
                        pairs.append((k, l))
                for l in range(nd, F):
                    m_ = self.odd[cl.nbr[l] - 2 * d]
                    j = cl.nbr[k] >> 1
                    wj = -1 if (m_ >> j) & 1 else 1
                    if (1 if cl.nbr[k] & 1 else -1) == wj:
                        pairs.append((k, l))
        else:
            pairs = [(k, l) for k in range(F) for l in range(k + 1, F)]
        cl.hinges = pairs
        self._cells[c] = cl
        return cl

    # ---- vertices (base units, rational)
    def vertices(self, c):
        v = self._verts.get(c)
        if v is not None:
            return v
        d, D = self.d, self.D
        if self.family == "simplex":
            out = []
            for j in range(D):
                if j != c:
                    x = [Fraction(0)] * D
                    x[j] = Fraction(1)
                    out.append(x)
        elif self.family == "orthoplex":
            out = []
            for i in range(d):
                x = [Fraction(0)] * d
                x[i] = Fraction(-1 if (c >> i) & 1 else 1)
                out.append(x)
        elif c < 2 * d:
            i = c >> 1
            a = 1 if c & 1 else -1
            out = []
            for m_ in range(2 ** d):
                if popcount(m_) % 2 == 0:
                    s_ = self.signs(m_)
                    if s_[i] == a:
                        out.append([Fraction(t, 2) for t in s_])
        else:
            w = self.signs(self.odd[c - 2 * d])
            out = []
            for k in range(d):
                x = [Fraction(t, 2) for t in w]
                x[k] = -x[k]
                out.append(x)
        self._verts[c] = out
        return out

    def centre_float(self, c):
        """The cell's centre at unit edge (floats), as the engine's container.ts computes it."""
        d, D = self.d, self.D
        o = [0.0] * D
        if self.family == "simplex":
            o = [1.0 / d] * D
            o[c] = 0.0
        elif self.family == "orthoplex":
            o = [(-1.0 if (c >> i) & 1 else 1.0) / d for i in range(d)]
        elif c < 2 * d:
            o[c >> 1] = 0.5 if c & 1 else -0.5
        else:
            w = self.signs(self.odd[c - 2 * d])
            o = [w[i] * (d - 2) / (2.0 * d) for i in range(d)]
        return [x / math.sqrt(2.0) for x in o]

    def normal_float(self, c):
        N, _, _ = self.normal(c)
        k = math.sqrt(self.kappa)
        return [to_float(x) / k for x in N]

    def frame_float(self, c):
        """The cell's frame: Gram-Schmidt twice on the projections of e_0, e_1, ... onto the tangent space, skipping residuals under 1e-9."""
        D = self.D
        cons = [self.normal_float(c)]
        if self.family == "simplex":
            cons.append([1.0 / math.sqrt(D)] * D)
        frame = []

        def proj(v, against):
            for b in against:
                t = sum(x * y for x, y in zip(b, v))
                for j in range(D):
                    v[j] -= t * b[j]

        for k in range(D):
            if len(frame) >= self.n:
                break
            v = [0.0] * D
            v[k] = 1.0
            for _ in range(2):
                proj(v, cons)
                proj(v, frame)
            nm = math.sqrt(sum(x * x for x in v))
            if nm < 1e-9:
                continue
            frame.append([x / nm for x in v])
        if len(frame) != self.n:
            raise NDError("frame failure")
        return frame

    def ambient_float(self, c, u, s):
        """s * centre + sum u_k frame_k in floats, as the engine's toAmbient."""
        o = self.centre_float(c)
        fr = self.frame_float(c)
        x = [s * t for t in o]
        for k in range(self.n):
            for j in range(self.D):
                x[j] += u[k] * fr[k][j]
        return x


_CONTS = {}


def get_cont(name):
    c = _CONTS.get(name)
    if c is None:
        c = _CONTS[name] = Cont(name)
    return c


# ------------------------------------------------------------------ developments, windows, joins


class Budget:
    def __init__(self, limit=None):
        self.steps = 0
        self.limit = MAX_STEPS if limit is None else limit

    def spend(self, k=1):
        self.steps += k
        if self.steps > self.limit:
            raise WorkLimit(f"the exact check needs more than {self.limit} steps")


class Node:
    __slots__ = ("cell", "parent", "via", "a")

    def __init__(self, cell, parent, via, a):
        self.cell, self.parent, self.via, self.a = cell, parent, via, a


def mirror(cn, cl, k):
    """The reflection vector m = N_c + N_g of the ridge between cell cl and its neighbour k, and m . m (cached)."""
    key = (cl.c, k)
    r = cn._mirrors.get(key)
    if r is None:
        m = [x + y for x, y in zip(cl.N, cl.Nn[k])]
        r = cn._mirrors[key] = (m, dot(m, m))
    return r


def reflect_dir(z, m, mm):
    c = 2 * dot(m, z) / mm
    return [zi - c * mi for zi, mi in zip(z, m)]


def develop(cn, c0, x, reach2, budget):
    """Half-chains from the point x (base units, in the hyperplane of cell c0) of reach sqrt(reach2): every chain a straight path of that length can follow."""
    nodes = [Node(c0, -1, -1, list(x))]
    kap = cn.kappa
    budget.spend()
    idx = 0
    while idx < len(nodes):
        nd = nodes[idx]
        cl = cn.cell(nd.cell)
        a = nd.a
        F = len(cl.nbr)
        parent_cell = nodes[nd.parent].cell if nd.parent >= 0 else -1
        S = [cl.hn[k] - dot(cl.Nn[k], a) for k in range(F)]
        for k in range(F):
            if cl.nbr[k] == parent_cell:
                continue
            Sk = S[k]
            if sg(Sk) <= 0:
                continue
            nu2 = cl.nu2[k]
            d2 = Sk * Sk / nu2
            if d2 > reach2:
                continue
            c1 = cl.Nc_dot[k]
            nu = [Nk - (c1 / kap) * Nc for Nk, Nc in zip(cl.Nn[k], cl.N)]
            f = Sk / nu2
            extra2 = 0
            for l in range(F):
                if l == k:
                    continue
                q = dot(cl.Nn[l], nu)
                viol = -S[l] + f * q
                if sg(viol) > 0:
                    tn2 = cl.nu2[l] - q * q / nu2
                    if sg(tn2) > 0:
                        e2 = viol * viol / tn2
                        if e2 > extra2:
                            extra2 = e2
            if d2 + extra2 > reach2:
                continue
            budget.spend()
            if len(nodes) >= MAX_NODES:
                raise WorkLimit(f"more than {MAX_NODES} chains from one centre")
            m, mm = mirror(cn, cl, k)
            q = cl.rp[k]
            ma = 2 * (dot(m, a) - dot(m, q)) / mm
            nodes.append(Node(cl.nbr[k], idx, k, [ai - ma * mi for ai, mi in zip(a, m)]))
        idx += 1
    return nodes


def in_window(cn, nodes, idx, w, beyond):
    """Is the direction w (in the frame of the end cell of node idx) in the closed window of the chain, and, if `beyond`, beyond its last ridge?"""
    z = list(w)
    cur = idx
    first = True
    while nodes[cur].parent >= 0:
        p = nodes[nodes[cur].parent]
        pc = cn.cell(p.cell)
        f = nodes[cur].via
        m, mm = mirror(cn, pc, f)
        z = reflect_dir(z, m, mm)
        F = len(pc.nbr)
        Sf = pc.hn[f] - dot(pc.Nn[f], p.a)
        nfz = dot(pc.Nn[f], z)
        if sg(nfz) < 0:
            return False
        if first and beyond and sg(nfz - Sf) < 0:
            return False
        first = False
        for k in range(F):
            if k == f:
                continue
            Sk = pc.hn[k] - dot(pc.Nn[k], p.a)
            if sg(Sf * dot(pc.Nn[k], z) - Sk * nfz) > 0:
                return False
        cur = nodes[cur].parent
    return True


def join_distance2(cn, A, B, limit2, budget, same, own=True):
    """The least squared length (base units) of a straight path joining the sources of the half-chain sets A and B, if <= limit2; else None.

    With `own` the sets are one ball's and a join of length exactly 0 is the trivial loop of a centre on a ridge, skipped. For two different
    balls (`own` False) length 0 is one point: an overlap, kept (4D.2)."""
    by_cell = {}
    for j, nd in enumerate(B):
        by_cell.setdefault(nd.cell, []).append(j)
    best = None
    for i, na in enumerate(A):
        for j in by_cell.get(na.cell, ()):
            if same and j <= i:
                continue
            w = [bj - ai for ai, bj in zip(na.a, B[j].a)]
            len2 = dot(w, w)
            if len2 > limit2 or (own and sg(len2) == 0) or (best is not None and len2 >= best):
                continue
            budget.spend()
            if not in_window(cn, A, i, w, True):
                continue
            if not in_window(cn, B, j, [-t for t in w], True):
                continue
            best = len2
    return best


# ------------------------------------------------------------------ the nearest hinge point inside a window


class Row:
    __slots__ = ("v", "b", "v2")

    def __init__(self, v, b):
        self.v = v
        self.b = b
        self.v2 = dot(v, v)


def window_rows(cn, nodes, idx):
    """The window of the chain as rows v . x <= b in the frame of the end cell, about the source a of the end node."""
    rows = []
    a_end = nodes[idx].a
    steps = []
    cur = idx
    while nodes[cur].parent >= 0:
        steps.append(cur)
        cur = nodes[cur].parent
    for i, cur in enumerate(steps):
        p = nodes[nodes[cur].parent]
        pc = cn.cell(p.cell)
        f = nodes[cur].via
        Sf = pc.hn[f] - dot(pc.Nn[f], p.a)

        def carry(v):
            for t in range(i, -1, -1):
                c2 = steps[t]
                q = nodes[nodes[c2].parent]
                m, mm = mirror(cn, cn.cell(q.cell), nodes[c2].via)
                v = reflect_dir(v, m, mm)
            return v

        v = carry([-t for t in pc.Nn[f]])
        rows.append(Row(v, dot(v, a_end)))
        for k in range(len(pc.nbr)):
            if k == f:
                continue
            Sk = pc.hn[k] - dot(pc.Nn[k], p.a)
            v = carry([Sf * x - Sk * y for x, y in zip(pc.Nn[k], pc.Nn[f])])
            rows.append(Row(v, dot(v, a_end)))
    return rows


def solve_linear(G, r):
    """Solves G x = r exactly (Gauss-Jordan); returns (x, G^-1) or None when G is singular."""
    m = len(G)
    A = [list(G[i]) + [Fraction(int(i == j)) for j in range(m)] for i in range(m)]
    for c in range(m):
        p = None
        for i in range(c, m):
            if sg(A[i][c]) != 0:
                p = i
                break
        if p is None:
            return None
        A[c], A[p] = A[p], A[c]
        pv = A[c][c]
        A[c] = [t / pv for t in A[c]]
        for i in range(m):
            if i != c and sg(A[i][c]) != 0:
                f = A[i][c]
                A[i] = [t - f * u for t, u in zip(A[i], A[c])]
    inv = [row[m:] for row in A]
    x = [sum((inv[i][j] * r[j] for j in range(m)), Fraction(0)) for i in range(m)]
    return x, inv


def nearest2(a, eqs, rows0, cap, budget):
    """The squared distance from a to { x : v_i . x = b_i (eqs), v_j . x <= b_j (rows) } if it is <= cap (else None); also the nearest point."""
    # a row the whole ball of squared radius cap about a satisfies cannot matter
    rows = []
    for r in rows0:
        margin = r.b - dot(r.v, a)
        if sg(margin) >= 0 and margin * margin >= cap * r.v2:
            continue
        rows.append(r)
    memo = {}
    active = []

    def solve(cap):
        key = tuple(sorted(active))
        if key in memo:
            return memo[key]
        budget.spend()
        E = list(eqs) + [rows[i] for i in active]
        m = len(E)
        G = [[dot(E[i].v, E[j].v) for j in range(m)] for i in range(m)]
        rr = [dot(E[i].v, a) - E[i].b for i in range(m)]
        sol = solve_linear(G, rr)
        if sol is None:
            memo[key] = None
            return None
        lam, inv = sol
        x = list(a)
        for i in range(m):
            if sg(lam[i]) != 0:
                x = [xj - lam[i] * vj for xj, vj in zip(x, E[i].v)]
        d2 = sum(((xj - aj) * (xj - aj) for xj, aj in zip(x, a)), Fraction(0))
        if d2 > cap:
            return None
        viol = []
        worst = 0
        for i, r in enumerate(rows):
            if i in active:
                continue
            v = dot(r.v, x) - r.b
            if sg(v) > 0:
                viol.append((v, i))
                ev = [dot(E[k].v, r.v) for k in range(m)]
                q = r.v2 - sum((ev[k] * inv[k][l] * ev[l] for k in range(m) for l in range(m)), Fraction(0))
                if sg(q) > 0:
                    t = v * v / q
                    if t > worst:
                        worst = t
        if d2 + worst > cap:
            return None
        out = None
        if not viol:
            out = (d2, x)
        else:
            viol.sort(key=lambda t: to_float(t[0]), reverse=True)
            for _, i in viol:
                active.append(i)
                r2 = solve(cap if out is None else min(cap, out[0]))
                active.pop()
                if r2 is not None and (out is None or r2[0] < out[0]):
                    out = r2
        if out is not None and out[0] <= cap:
            memo[key] = out
        return out

    return solve(cap)


def hinge_dist2(cn, nodes, idx, reach2, budget):
    """The least squared distance from the source to a hinge point of the cell of node idx seen through its window (<= reach2), or None."""
    nd = nodes[idx]
    cl = cn.cell(nd.cell)
    a = nd.a
    F = len(cl.nbr)
    S = [cl.hn[k] - dot(cl.Nn[k], a) for k in range(F)]
    wrows = None
    best = None
    kap = cn.kappa
    for (k, l) in cl.hinges:
        if S[k] * S[k] > reach2 * cl.nu2[k] or S[l] * S[l] > reach2 * cl.nu2[l]:
            continue
        # a lower bound: the distance to the hinge's flat inside the hyperplane
        c = dot(cl.Nn[k], cl.Nn[l]) - cl.Nc_dot[k] * cl.Nc_dot[l] / kap
        det = cl.nu2[k] * cl.nu2[l] - c * c
        if sg(det) <= 0:
            continue
        lb = (S[k] * S[k] * cl.nu2[l] - 2 * S[k] * S[l] * c + S[l] * S[l] * cl.nu2[k]) / det
        if lb > reach2:
            continue
        budget.spend()
        if wrows is None:
            wrows = window_rows(cn, nodes, idx)
        eqs = [Row(cl.N, cl.h)]
        if cn.family == "simplex":
            eqs.append(Row([Fraction(1)] * cn.D, Fraction(1)))
        eqs.append(Row(cl.Nn[k], cl.hn[k]))
        eqs.append(Row(cl.Nn[l], cl.hn[l]))
        rows = [Row(cl.Nn[j], cl.hn[j]) for j in range(F) if j != k and j != l] + wrows
        cap = reach2 if best is None else min(reach2, best[0])
        r = nearest2(a, eqs, rows, cap, budget)
        if r is not None and (best is None or r[0] < best[0]):
            best = r
    return None if best is None else best[0]


# ------------------------------------------------------------------ configurations


SQRT2_Q = Fraction(sqrt_floor(2 * 10 ** 120), 10 ** 60)  # sqrt 2 to 60 digits


def fmt(x):
    return f"{x:.10g}"


class Cfg:
    pass


def parse_config(polytope, P, config, to_fraction):
    cn = get_cont(polytope)
    if not isinstance(config, dict):
        raise NDError("the configuration must be an object with s and pieces")
    s = to_fraction(config.get("s"), "s")
    if not 0 < s <= MAX_S:
        raise NDError("s must be greater than 0 and at most 1e4")
    pieces = config.get("pieces")
    if not isinstance(pieces, list):
        raise NDError("the configuration needs a list of pieces")
    if len(pieces) != P:
        raise NDError(f"the configuration has {len(pieces)} pieces; the problem needs {P}")
    cells, U = [], []
    for n_, p in enumerate(pieces):
        if not isinstance(p, dict):
            raise NDError(f"piece {n_} must be an object with cell and u")
        extra = set(p) - {"cell", "u", "q"}
        if extra:
            raise NDError(f"piece {n_} has unknown fields {sorted(extra)}")
        if "q" in p:
            raise NDError(f"piece {n_} is a ball and must not have a rotation q")
        c = p.get("cell")
        if isinstance(c, bool) or not isinstance(c, int) or not 0 <= c < cn.C:
            raise NDError(f"piece {n_}: cell must be a whole number from 0 to {cn.C - 1}")
        u = p.get("u")
        if not isinstance(u, list) or len(u) != cn.n:
            raise NDError(f"piece {n_}: u must be {cn.n} numbers")
        uu = [to_fraction(x, f"piece {n_} u[{k}]") for k, x in enumerate(u)]
        if any(abs(x) > MAX_S for x in uu):
            raise NDError(f"piece {n_}: u is out of range")
        cells.append(c)
        U.append(uu)
    cfg = Cfg()
    cfg.cont, cfg.polytope, cfg.s, cfg.P, cfg.cells, cfg.U = cn, polytope, s, P, cells, U
    return cfg


def independent_vertices(cn, c):
    """n + 1 affinely independent vertices of the cell (exact), greedy in the list's order."""
    vs = cn.vertices(c)
    base = vs[0]
    chosen = [base]
    basis = []  # reduced row vectors of differences

    def reduce(v):
        v = list(v)
        for b, piv in basis:
            if v[piv] != 0:
                f = v[piv] / b[piv]
                v = [x - f * y for x, y in zip(v, b)]
        return v

    for v in vs[1:]:
        r = reduce([x - y for x, y in zip(v, base)])
        piv = next((i for i, t in enumerate(r) if t != 0), None)
        if piv is None:
            continue
        basis.append((r, piv))
        chosen.append(v)
        if len(chosen) == cn.n + 1:
            break
    if len(chosen) != cn.n + 1:
        raise NDError("cell vertices do not span the cell")
    return chosen


def read_centre(cn, c, u, s):
    """The exact base-unit point for the engine's (cell c, u at edge s): an affine combination of n + 1 vertices of the cell with dyadic weights."""
    xf = cn.ambient_float(c, [float(t) if not isinstance(t, float) else t for t in u], float(s))
    y = [Fraction(t) * SQRT2_Q / s for t in xf]
    vs = independent_vertices(cn, c)
    v0 = vs[0]
    B = [[a - b for a, b in zip(v, v0)] for v in vs[1:]]
    n = cn.n
    rhs = [x - b for x, b in zip(y, v0)]
    G = [[dot(B[i], B[j]) for j in range(n)] for i in range(n)]
    r = [dot(B[i], rhs) for i in range(n)]
    sol = solve_linear(G, r)
    if sol is None:
        raise NDError("singular reading")
    w = sol[0]
    q = Fraction(1, 2 ** 60)
    wd = [Fraction(round(x / q)) * q for x in w]
    x = list(v0)
    for k in range(n):
        x = [xi + wd[k] * bi for xi, bi in zip(x, B[k])]
    return x


def reading_error(polytope, c, u, s):
    """For the documentation's measurement: the distance (relative to s) between the engine's float point for (cell c, u at edge s) and the exact
    point read_centre gives (before any clamping)."""
    cn = get_cont(polytope)
    s = Fraction(s)
    xf = cn.ambient_float(c, [float(t) for t in u], float(s))
    x = read_centre(cn, c, [Fraction(t) for t in u], s)
    err2 = sum(((Fraction(a) - b * s / SQRT2_Q) ** 2 for a, b in zip(xf, x)), Fraction(0))
    return math.sqrt(float(err2)) / float(s)


def clamp_exact(cn, c, x, s):
    """Moves x onto the facet planes of cell c that it violates, if by less than SNAP s; returns the point and the real distance moved (squared, base units)."""
    cl = cn.cell(c)
    kap = cn.kappa
    x0 = list(x)
    lim2 = (SNAP * s) ** 2 * 2 / (s * s)  # (SNAP s)^2 in base units^2: base = real sqrt2 / s
    for _ in range(8):
        moved = False
        for k in range(len(cl.nbr)):
            S = cl.hn[k] - dot(cl.Nn[k], x)
            if sg(S) < 0:
                nu = [Nk - (cl.Nc_dot[k] / kap) * Nc for Nk, Nc in zip(cl.Nn[k], cl.N)]
                f = S / cl.nu2[k]
                x = [xi + f * ni for xi, ni in zip(x, nu)]
                moved = True
        if not moved:
            break
    mv = dot([a - b for a, b in zip(x, x0)], [a - b for a, b in zip(x, x0)])
    return x, mv


# ------------------------------------------------------------------ verdicts


class Verdict:
    def __init__(self, cfg, strict):
        self.cfg = cfg
        self.strict = strict
        self.reasons = []
        self.margins = {}
        self.stats = {}

    @property
    def ok(self):
        return not self.reasons

    def fail(self, msg):
        if len(self.reasons) < 50:
            self.reasons.append(msg)


def sqrt2_lo_pow(n):
    """(sqrt 2)^n rounded down."""
    r = Fraction(sqrt_floor(2 * 10 ** 120), 10 ** 60)
    return Fraction(2 ** (n // 2)) * (r if n % 2 else 1)


def sqrt_hi(x):
    """sqrt(x) rounded up (a rational), x a positive integer."""
    return Fraction(sqrt_floor(x * 10 ** 120) + 1, 10 ** 60)


def boundary_volume_hi(cn):
    """A rational upper bound of the n-volume of the boundary at unit edge: a regular n-simplex has volume sqrt(n+1) / (n! 2^(n/2)), an
    n-demicube 2^(-n/2) (1 - 2^(n-1)/n!); the boundary has (d+1), 2^d simplices, or 2d demicubes and 2^(d-1) simplices."""
    d, n = cn.d, cn.n
    nf = math.factorial(n)
    r2n = sqrt2_lo_pow(n)
    simplex = sqrt_hi(d) / (nf * r2n)
    if cn.family == "simplex":
        return (d + 1) * simplex
    if cn.family == "orthoplex":
        return 2 ** d * simplex
    return 2 * d * (1 - Fraction(2 ** (n - 1), nf)) / r2n + 2 ** (d - 1) * simplex


def ball_volume_lo(n):
    """A rational lower bound of the volume of a ball of radius 1/2 in n dimensions."""
    num, den, k = {4: (1, 2, 2), 5: (8, 15, 2), 6: (1, 6, 3), 7: (16, 105, 3), 8: (1, 24, 4), 9: (32, 945, 4)}[n]
    return Fraction(num, den) * PI_LO ** k / 2 ** n


def ball_volume_ok(cn, P, s):
    """P balls fit in the boundary's volume, with pi and the square roots rounded the safe way: False is certain."""
    return P * ball_volume_lo(cn.n) <= boundary_volume_hi(cn) * s ** cn.n


def check(polytope, piece, P, config, to_fraction, strict=True, extra=0):
    """Decides validity exactly. `extra` (a length, site units) widens the search past the thresholds so that margins up to that much are found
    and reported (the fixtures' builder asks for 2e-12 to show their margins are at least 1e-12); the verdict does not depend on it."""
    if polytope not in NAMES:
        raise NDError(f"unknown container {polytope!r}")
    if piece != "ball":
        raise NDError(f"the {polytope} takes balls only, not {piece!r} (BOUNDARY.md section 4D)")
    if isinstance(P, bool) or not isinstance(P, int) or not 1 <= P <= MAX_P:
        raise NDError(f"P must be a whole number from 1 to {MAX_P}")
    cfg = parse_config(polytope, P, config, to_fraction)
    v = Verdict(cfg, strict)
    cn = cfg.cont
    s = cfg.s
    if not ball_volume_ok(cn, P, s):
        v.fail(f"{P} balls have more volume than the boundary of the {polytope} at s = {fmt(float(s))}")
        return v
    try:
        _check_balls(cn, cfg, v, Fraction(extra))
    except WorkLimit as e:
        raise NDError(str(e)) from None
    return v


def _check_balls(cn, cfg, v, extra=Fraction(0)):
    s, strict = cfg.s, v.strict
    scale2 = Fraction(2) / (s * s)  # real^2 -> base^2
    sp = READING_PAIR * s if strict else Fraction(0)
    sc = READING_CLEAR * s if strict else Fraction(0)
    pair_lim2 = (1 + sp) ** 2 * scale2
    clr_lim2 = (Fraction(1, 2) + sc) ** 2 * scale2
    half = max(Fraction(1, 2) + sc, (1 + sp) / 2) + extra
    reach2 = half * half * scale2
    # pairs and loops are searched up to a length of 2 * half, so that margins up to 2 extra show
    pair_search2 = (2 * half) ** 2 * scale2
    budget = Budget()
    pts = []
    for n_, (c, u) in enumerate(zip(cfg.cells, cfg.U)):
        x = read_centre(cn, c, u, s)
        x, mv = clamp_exact(cn, c, x, s)
        cl = cn.cell(c)
        if any(sg(cl.hn[k] - dot(cl.Nn[k], x)) < 0 for k in range(len(cl.nbr))):
            v.fail(f"piece {n_}: its centre lies outside its anchor cell {c}")
            continue
        if mv * s * s / 2 > (SNAP * s) ** 2:
            v.fail(f"piece {n_}: its centre lies outside its anchor cell {c} by more than 1e-12 of an edge")
            continue
        anchors = [(c, x)]
        for k in range(len(cl.nbr)):
            if sg(cl.hn[k] - dot(cl.Nn[k], x)) == 0:
                # on the ridge: it belongs to the neighbour too (x lies in its hyperplane exactly)
                anchors.append((cl.nbr[k], x))
        pts.append((n_, x, anchors))
    if not v.ok:
        return
    # clearance: a centre nearer than the threshold to a hinge of its own cell is refused before any chain is built
    chains = {}
    clear = []
    for n_, x, anchors in pts:
        best = None
        for (c, xa) in anchors:
            root = develop(cn, c, xa, Fraction(0), budget)
            h2 = hinge_dist2(cn, root, 0, reach2, budget)
            if h2 is not None and (best is None or h2 < best):
                best = h2
        bad = best is not None and (best < clr_lim2 or (strict and best == clr_lim2))
        if not bad:
            cs = []
            for (c, xa) in anchors:
                nodes = develop(cn, c, xa, reach2, budget)
                cs.append(nodes)
                for i in range(len(nodes)):
                    h2 = hinge_dist2(cn, nodes, i, reach2, budget)
                    if h2 is not None and (best is None or h2 < best):
                        best = h2
            chains[n_] = cs
            bad = best is not None and (best < clr_lim2 or (strict and best == clr_lim2))
        clear.append(best)
        if bad:
            v.fail(f"ball {n_} is {fmt(math.sqrt(to_float(best / scale2)))} from a (d-3)-face of the {cfg.polytope}; it must be {'more than' if strict else 'at least'} 1/2")
    mins = [b for b in clear if b is not None]
    v.margins["clearance_sq"] = None if not mins else min(mins) / scale2
    if not v.ok:
        return
    # pairs (the chord shortcut: R^D distance)
    best_pair = None
    npairs = 0
    for i in range(len(pts)):
        for j in range(i + 1, len(pts)):
            xi, xj = pts[i][1], pts[j][1]
            diff = [p - q for p, q in zip(xi, xj)]
            if dot(diff, diff) > pair_search2:
                continue
            npairs += 1
            m = None
            if sg(dot(diff, diff)) == 0:
                # one point: two different balls at the same place overlap whatever the chains say (4D.2)
                m = Fraction(0)
            for A in chains[pts[i][0]]:
                for B in chains[pts[j][0]]:
                    r = join_distance2(cn, A, B, pair_search2, budget, False, False)
                    if r is not None and (m is None or r < m):
                        m = r
            if m is None:
                continue
            best_pair = m if best_pair is None or m < best_pair else best_pair
            if m < pair_lim2 or (strict and m == pair_lim2):
                v.fail(f"balls {pts[i][0]} and {pts[j][0]} overlap: a path of length {fmt(math.sqrt(to_float(m / scale2)))} joins their centres; they must be {'more than' if strict else 'at least'} 1 apart")
    best_loop = None
    for n_, x, anchors in pts:
        cs = chains[n_]
        m = None
        for a in range(len(cs)):
            for b in range(a, len(cs)):
                r = join_distance2(cn, cs[a], cs[b], pair_search2, budget, a == b)
                if r is not None and (m is None or r < m):
                    m = r
        if m is None:
            continue
        best_loop = m if best_loop is None or m < best_loop else best_loop
        if m < pair_lim2 or (strict and m == pair_lim2):
            v.fail(f"ball {n_} overlaps itself: a geodesic loop of length {fmt(math.sqrt(to_float(m / scale2)))} leaves and returns to its centre")
    v.margins["distance_sq"] = None if best_pair is None else best_pair / scale2
    v.margins["loop_sq"] = None if best_loop is None else best_loop / scale2
    v.stats["pairs_within_1_in_Rd"] = npairs
    v.stats["chains"] = sum(len(nodes) for cs in chains.values() for nodes in cs)
    v.stats["steps"] = budget.steps


def read_pieces(cn, cfg):
    """The exact centres of a configuration (clamped into their closed cells) with the cells each one belongs to; raises NDError if one lies outside."""
    s = cfg.s
    pts = []
    for n_, (c, u) in enumerate(zip(cfg.cells, cfg.U)):
        x = read_centre(cn, c, u, s)
        x, mv = clamp_exact(cn, c, x, s)
        cl = cn.cell(c)
        if any(sg(cl.hn[k] - dot(cl.Nn[k], x)) < 0 for k in range(len(cl.nbr))) or mv * s * s / 2 > (SNAP * s) ** 2:
            raise NDError(f"piece {n_}: its centre lies outside its anchor cell {c}")
        anchors = [(c, x)] + [(cl.nbr[k], x) for k in range(len(cl.nbr)) if sg(cl.hn[k] - dot(cl.Nn[k], x)) == 0]
        pts.append((n_, x, anchors))
    return pts


def measure(polytope, P, config, to_fraction, reach=Fraction(1, 2) * (1 + Fraction(1, 10 ** 12))):
    """For cross-checks: the exact squared distances (site units) the rule is made of, without a verdict: per centre the geodesic distance to the
    hinges when it is at most `reach` (else None), per pair of centres closer than 1 in R^D the shortest straight path when it is below 1
    (else None), per centre the shortest loop below 1. A centre within 1/2 of a hinge of its own cell is not developed further (pairs: the chord)."""
    cfg = parse_config(polytope, P, config, to_fraction)
    cn, s = cfg.cont, cfg.s
    scale2 = Fraction(2) / (s * s)
    reach2 = reach * reach * scale2
    pair_lim2 = scale2
    budget = Budget()
    pts = read_pieces(cn, cfg)
    clear, chains, crowded = [], {}, {}
    for n_, x, anchors in pts:
        best = None
        for (c, xa) in anchors:
            nodes = develop(cn, c, xa, Fraction(0), budget)
            h2 = hinge_dist2(cn, nodes, 0, reach2, budget)
            if h2 is not None and (best is None or h2 < best):
                best = h2
        crowded[n_] = best is not None and best < Fraction(1, 4) * scale2
        if not crowded[n_]:
            cs = []
            for (c, xa) in anchors:
                nodes = develop(cn, c, xa, reach2, budget)
                cs.append(nodes)
                for i in range(len(nodes)):
                    h2 = hinge_dist2(cn, nodes, i, reach2, budget)
                    if h2 is not None and (best is None or h2 < best):
                        best = h2
            chains[n_] = cs
        clear.append(None if best is None else best / scale2)
    pairs = {}
    for i in range(len(pts)):
        for j in range(i + 1, len(pts)):
            diff = [p - q for p, q in zip(pts[i][1], pts[j][1])]
            q2 = dot(diff, diff)
            if q2 >= pair_lim2:
                continue
            if crowded[i] or crowded[j]:
                pairs[(i, j)] = q2 / scale2
                continue
            m = None
            for A in chains[i]:
                for B in chains[j]:
                    r = join_distance2(cn, A, B, pair_lim2 - Fraction(1, 10 ** 30), budget, False, False)
                    if r is not None and (m is None or r < m):
                        m = r
            pairs[(i, j)] = None if m is None else m / scale2
    loops = []
    for n_, x, anchors in pts:
        if crowded[n_]:
            loops.append(None)
            continue
        cs = chains[n_]
        m = None
        for a in range(len(cs)):
            for b in range(a, len(cs)):
                r = join_distance2(cn, cs[a], cs[b], pair_lim2 - Fraction(1, 10 ** 30), budget, a == b)
                if r is not None and (m is None or r < m):
                    m = r
        loops.append(None if m is None else m / scale2)
    return {"clearance_sq": clear, "pairs_sq": pairs, "loops_sq": loops, "steps": budget.steps}


def sqrt_down(x):
    """The largest binary64 whose square is <= x >= 0 (a rational, or a number of Q(sqrt k): its lower bound)."""
    if hasattr(x, "bounds"):
        x = x.bounds()[0]
    if x <= 0:
        return 0.0
    f = math.sqrt(float(x))
    while Fraction(f) ** 2 > x:
        f = math.nextafter(f, 0.0)
    while Fraction(math.nextafter(f, math.inf)) ** 2 <= x:
        f = math.nextafter(f, math.inf)
    return f


def summary(v):
    """Readable margins of a verdict (site lengths, rounded down so that they are never overstated)."""
    m = v.margins
    out = {"checker": CHECKER, "version": VERSION, "strict": v.strict, "valid": v.ok}
    out["min_clearance"] = None if m.get("clearance_sq") is None else sqrt_down(m["clearance_sq"])
    out["min_distance"] = None if m.get("distance_sq") is None else sqrt_down(m["distance_sq"])
    out["min_loop"] = None if m.get("loop_sq") is None else sqrt_down(m["loop_sq"])
    out["chains"] = v.stats.get("chains")
    out["steps"] = v.stats.get("steps")
    return out


def exact_bounds(polytope, P):
    """The proven lower bounds on s (BOUNDARY.md 4D.5) as (kind, text, check), check(x) deciding exactly that the rational x is not above the bound."""
    cn = get_cont(polytope)
    d, n = cn.d, cn.n
    ball_lo = ball_volume_lo(n)
    vol_hi = boundary_volume_hi(cn)

    def volume_check(x):
        return vol_hi * x ** n <= P * ball_lo

    def one_check(x):
        if cn.family == "demicube":
            return x * x * (n * n - 5 * n + 8) <= 2 * (n - 1)
        return 4 * x * x <= d * (d - 2)

    return [("volume", "(P V_(d-1) 2^-(d-1) / vol_1)^(1/(d-1))", volume_check), ("one piece", "1 / (2 rho)", one_check)]


def main(argv):
    import json

    if len(argv) != 2:
        print("usage: python3 python/boundary/ndpoly.py boundary/<d><family>/ball/<P> <config.json>")
        return 2
    m = ID_RE.match(argv[0])
    if not m:
        print(f"FAIL {argv[0]!r} is not an nD container problem id")
        return 1
    with open(argv[1]) as fh:
        config = json.load(fh)
    sys.path.insert(0, __file__.rsplit("/", 1)[0])
    import certify

    try:
        print(json.dumps(certify.certify(argv[0], config), indent=2))
    except certify.CertError as e:
        print(f"FAIL {e}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
