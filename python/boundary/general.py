"""Exact certifier for balls on the boundary of a regular 4-polytope ("3D on 4D", the 5-cell, 16-cell, 24-cell, 120-cell,
600-cell and also the tesseract), with developments along chains of cells.

Problem boundary/<polytope>/ball/<P>: P balls of diameter 1 on the boundary of the polytope of edge s; the score is s.
The conventions are those of packages/boundary/src/types.ts and docs/spaces/BOUNDARY.md: a configuration is
{s, pieces: [{cell, u}]}, cell an index into the polytope's cells in the engine's order (facets sorted by the outward
unit normal with compareNormals), u the centre in that cell's local frame at the configuration's s (the engine's frame:
Gram-Schmidt of the standard basis projected on the cell's hyperplane), the cell's centre at the origin.

Exactness. Everything below is decided in exact rational arithmetic, or in Q(sqrt 5) for the 120-cell and the
600-cell, with integer square roots only to report numbers. The polytope is built from exact vertices, a copy of the
unit-edge polytope up to a scale (BASE: edge e0, e0^2 exact; the 5-cell lives in R^5 as the simplex with vertices
e_i - (1/5) 1, edge^2 = 2, the 16-cell and 24-cell in Z^4 with edge^2 = 2, the 600-cell with circumradius 1 and the
120-cell with edge 3 - sqrt 5, both in Q(sqrt 5)^4). Lengths of the base polytope convert to the problem's by the
factor s / e0 (a ball of diameter 1 has diameter e0 / s in base units, so every test compares a squared base length
with e0^2 / s^2, an exact number of the field).

The centre of a ball is read exactly from the engine's numbers: the ambient point s * origin + sum u_k frame_k of the
cell (floats, as the engine computes it) is written as an affine combination of four independent vertices of the
cell with dyadic weights (rounded to 2^-60 and summing to 1), and the exact point is the same combination of the
exact vertices. That point lies in the cell's hyperplane exactly and differs from the engine's by about 1e-16 s; a
point that falls outside its closed cell by less than 1e-12 s is moved onto the facet planes it violates (as the
engine's clampToCell does) before the exact test "the centre lies in its closed cell".

Charts. A cell's hyperplane has an affine chart: lambda in the field^3 with x = O + lambda B, O a vertex of the cell
and B three of its edge vectors (rational or in the field). Distances need the metric G = B B^T (3 x 3, exact); the
windows, gates and orderings use only affine operations (cross products of chart vectors are covariant normals). The
gluing across a ridge R between cells c and c' is the rotation of R^n in the plane of the two outward normals that
fixes R; with N, N' the normals (equal norms a, k = N.N'/a) it is the exact formula
    Rot x = x + (k - 1) [ (N.x) N / a + (E.x) E / ((1 - k^2) a) ] + [ (N.x) E - (E.x) N ] / a,   E = N' - k N,
which takes N to N' and fixes the orthogonal complement of N, N'; no square root appears. The development of a cell
along a chain is an affine map of charts, a composition of the inverse gluings.

Chains and windows (BOUNDARY.md section 4A). From a centre u of a cell the straight paths that avoid the hinges (the
edges) cross a chain of ridges; the set of directions that follow a chain is the cone from u over the gate (the last
ridge clipped by the previous cone). The enumerator computes these cones exactly (Sutherland-Hodgman in the field),
and keeps a chain when (1) u lies strictly on the inner side of the ridge's plane, (2) the clipped gate has positive
area, (3) the distance from u to the gate is at most the reach, here the ball's diameter. A ball is valid exactly when
its clearance, every pair distance and every loop are at least the right numbers (4.2); each is the least over the
chains, tested on the exact image of the other centre, which must lie in the closed cone.

Usage: imported by certify.py. python3 python/boundary/general.py boundary/<polytope>/ball/<P> <config.json> certifies one configuration
(strictly valid, or the reasons) and prints its margins.
"""

import itertools
import math
from fractions import Fraction

# ------------------------------------------------------------------ the field Q(sqrt 5)


def sgn(x):
    """-1, 0 or 1: the exact sign of a Fraction, an int or a Q5."""
    if isinstance(x, Q5):
        return x.sign()
    return (x > 0) - (x < 0)


class Q5:
    """(n1 + n2 sqrt 5) / d with integers n1, n2 and d > 0 in lowest terms (gcd(n1, n2, d) = 1), so equal numbers have equal
    fields. Integer arithmetic with one gcd per operation is several times faster than pairs of Fractions."""

    __slots__ = ("n1", "n2", "d")

    def __init__(self, a=0, b=0, _raw=None):
        if _raw is not None:
            self.n1, self.n2, self.d = _raw
            return
        a = a if isinstance(a, Fraction) else Fraction(a)
        b = b if isinstance(b, Fraction) else Fraction(b)
        d = a.denominator * b.denominator // math.gcd(a.denominator, b.denominator)
        n1, n2 = a.numerator * (d // a.denominator), b.numerator * (d // b.denominator)
        g = math.gcd(n1, n2, d)
        self.n1, self.n2, self.d = n1 // g, n2 // g, d // g

    @property
    def a(self):
        return Fraction(self.n1, self.d)

    @property
    def b(self):
        return Fraction(self.n2, self.d)

    @staticmethod
    def _mk(n1, n2, d):
        g = math.gcd(n1, n2, d)
        if g != 1:
            n1, n2, d = n1 // g, n2 // g, d // g
        return Q5(_raw=(n1, n2, d))

    def __add__(self, o):
        if isinstance(o, Q5):
            return Q5._mk(self.n1 * o.d + o.n1 * self.d, self.n2 * o.d + o.n2 * self.d, self.d * o.d)
        if isinstance(o, int):
            return Q5._mk(self.n1 + o * self.d, self.n2, self.d)
        return Q5._mk(self.n1 * o.denominator + o.numerator * self.d, self.n2 * o.denominator, self.d * o.denominator)

    __radd__ = __add__

    def __neg__(self):
        return Q5(_raw=(-self.n1, -self.n2, self.d))

    def __sub__(self, o):
        if isinstance(o, Q5):
            return Q5._mk(self.n1 * o.d - o.n1 * self.d, self.n2 * o.d - o.n2 * self.d, self.d * o.d)
        if isinstance(o, int):
            return Q5._mk(self.n1 - o * self.d, self.n2, self.d)
        return Q5._mk(self.n1 * o.denominator - o.numerator * self.d, self.n2 * o.denominator, self.d * o.denominator)

    def __rsub__(self, o):
        return (-self) + o

    def __mul__(self, o):
        if isinstance(o, Q5):
            return Q5._mk(self.n1 * o.n1 + 5 * self.n2 * o.n2, self.n1 * o.n2 + self.n2 * o.n1, self.d * o.d)
        if isinstance(o, int):
            return Q5._mk(self.n1 * o, self.n2 * o, self.d)
        return Q5._mk(self.n1 * o.numerator, self.n2 * o.numerator, self.d * o.denominator)

    __rmul__ = __mul__

    def inv(self):
        n = self.n1 * self.n1 - 5 * self.n2 * self.n2
        if n == 0:
            raise ZeroDivisionError("Q5 inverse of 0")
        # d / (n1 + n2 s) = d (n1 - n2 s) / n
        if n < 0:
            return Q5._mk(-self.d * self.n1, self.d * self.n2, -n)
        return Q5._mk(self.d * self.n1, -self.d * self.n2, n)

    def __truediv__(self, o):
        if isinstance(o, Q5):
            return self * o.inv()
        if isinstance(o, int):
            o = Fraction(o)
        return Q5._mk(self.n1 * o.denominator, self.n2 * o.denominator, self.d * o.numerator) if o.numerator > 0 else Q5._mk(-self.n1 * o.denominator, -self.n2 * o.denominator, -self.d * o.numerator)

    def __rtruediv__(self, o):
        return self.inv() * o

    def sign(self):
        sa, sb = (self.n1 > 0) - (self.n1 < 0), (self.n2 > 0) - (self.n2 < 0)
        if sb == 0:
            return sa
        if sa == 0 or sa == sb:
            return sb
        # opposite signs: compare n1^2 with 5 n2^2 (never equal for nonzero integers)
        return sa if self.n1 * self.n1 > 5 * self.n2 * self.n2 else sb

    def __eq__(self, o):
        if isinstance(o, Q5):
            return self.n1 == o.n1 and self.n2 == o.n2 and self.d == o.d
        if isinstance(o, int):
            return self.n2 == 0 and self.n1 == o * self.d
        return self.n2 == 0 and self.n1 * o.denominator == o.numerator * self.d

    def __ne__(self, o):
        return not self == o

    def __hash__(self):
        return hash((self.n1, self.n2, self.d)) if self.n2 else hash(Fraction(self.n1, self.d))

    def __lt__(self, o):
        return (self - o).sign() < 0

    def __le__(self, o):
        return (self - o).sign() <= 0

    def __gt__(self, o):
        return (self - o).sign() > 0

    def __ge__(self, o):
        return (self - o).sign() >= 0

    def __repr__(self):
        return f"(({self.n1} + {self.n2}*sqrt5)/{self.d})"

    def bounds(self):
        """Rational (lo, hi) with lo <= value <= hi, 1e-40 apart."""
        lo, hi = SQRT5_LO, SQRT5_HI
        a, b = Fraction(self.n1, self.d), Fraction(self.n2, self.d)
        if b >= 0:
            return a + b * lo, a + b * hi
        return a + b * hi, a + b * lo

    def __float__(self):
        if self.n2 == 0:
            return self.n1 / self.d
        lo, hi = self.bounds()
        return float((lo + hi) / 2)


_S = 10 ** 45
SQRT5_LO = Fraction(math.isqrt(5 * _S * _S), _S)
SQRT5_HI = SQRT5_LO + Fraction(1, _S)


def lower_float(x):
    """The largest binary64 not above the exact (rational or Q5) x."""
    lo = x.bounds()[0] if isinstance(x, Q5) else Fraction(x)
    f = float(lo)
    while Fraction(f) > lo:
        f = math.nextafter(f, -math.inf)
    return f


def upper_float(x):
    hi = x.bounds()[1] if isinstance(x, Q5) else Fraction(x)
    f = float(hi)
    while Fraction(f) < hi:
        f = math.nextafter(f, math.inf)
    return f


# ------------------------------------------------------------------ small linear algebra over a field (lists)


def vdot(a, b):
    s = a[0] * b[0]
    for i in range(1, len(a)):
        s = s + a[i] * b[i]
    return s


def vsub(a, b):
    return [x - y for x, y in zip(a, b)]


def vadd(a, b):
    return [x + y for x, y in zip(a, b)]


def vscale(k, a):
    return [k * x for x in a]


def cross(a, b):
    return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]


def det3(m):
    return (
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    )


def inv3(m):
    """Inverse of a 3x3 matrix over a field (adjugate / determinant)."""
    d = det3(m)
    if sgn(d) == 0:
        raise ValueError("singular")
    c = [
        [m[1][1] * m[2][2] - m[1][2] * m[2][1], m[0][2] * m[2][1] - m[0][1] * m[2][2], m[0][1] * m[1][2] - m[0][2] * m[1][1]],
        [m[1][2] * m[2][0] - m[1][0] * m[2][2], m[0][0] * m[2][2] - m[0][2] * m[2][0], m[0][2] * m[1][0] - m[0][0] * m[1][2]],
        [m[1][0] * m[2][1] - m[1][1] * m[2][0], m[0][1] * m[2][0] - m[0][0] * m[2][1], m[0][0] * m[1][1] - m[0][1] * m[1][0]],
    ]
    di = 1 / d if not isinstance(d, Q5) else d.inv()
    return [[c[i][j] * di for j in range(3)] for i in range(3)]


def matvec(m, v):
    return [vdot(row, v) for row in m]


def matmat(a, b):
    cols = [[b[i][j] for i in range(3)] for j in range(3)]
    out = [[None] * 3 for _ in range(3)]
    for j, col in enumerate(cols):
        for i in range(3):
            out[i][j] = vdot(a[i], col)
    return out


# ------------------------------------------------------------------ the polytopes' exact vertices

NAMES = ("5cell", "16cell", "24cell", "120cell", "600cell", "tesseract")

EVEN_PERMUTATIONS = [p for p in itertools.permutations(range(4)) if sum(1 for i in range(4) for j in range(i + 1, 4) if p[i] > p[j]) % 2 == 0]
ALL_PERMUTATIONS = list(itertools.permutations(range(4)))


def signed_points(base, perms, zero):
    """Python twin of containers.ts signedPoints: the points point[perm[i]] = base[i], with every sign choice on the
    nonzero entries (the first nonzero entry's sign changing slowest, + before -), in order of first appearance."""
    nz = [i for i, x in enumerate(base) if sgn(x) != 0]
    out, seen = [], set()
    for perm in perms:
        for m in range(1 << len(nz)):
            p = [zero] * 4
            for i, x in enumerate(base):
                k = nz.index(i) if i in nz else -1
                neg = k >= 0 and (m >> (len(nz) - 1 - k)) & 1
                p[perm[i]] = -x if neg else x
            key = tuple(p)
            if key not in seen:
                seen.add(key)
                out.append(p)
    return out


def raw_vertices(name):
    """(vertices, K) of the base polytope in the order of containers.ts rawVertices; K is 'Q' or 'Q5'."""
    F = Fraction
    if name == "5cell":
        c = F(1, 5)
        return [[(F(1) if i == j else F(0)) - c for j in range(5)] for i in range(5)], "Q"
    if name == "16cell":
        out = []
        for i in range(4):
            for sg in (-1, 1):
                p = [F(0)] * 4
                p[i] = F(sg)
                out.append(p)
        return out, "Q"
    if name == "24cell":
        return signed_points([F(1), F(1), F(0), F(0)], ALL_PERMUTATIONS, F(0)), "Q"
    if name == "tesseract":
        return [[F(1, 2) if (v >> k) & 1 else F(-1, 2) for k in range(4)] for v in range(16)], "Q"
    zero = Q5(0)
    one = Q5(1)
    sq5 = Q5(0, 1)
    phi = Q5(F(1, 2), F(1, 2))
    iphi = Q5(F(-1, 2), F(1, 2))  # 1 / phi
    half = Q5(F(1, 2))
    if name == "600cell":
        v = signed_points([one, zero, zero, zero], ALL_PERMUTATIONS, zero)
        v += signed_points([half, half, half, half], [(0, 1, 2, 3)], zero)
        v += signed_points([phi * half, half, iphi * half, zero], EVEN_PERMUTATIONS, zero)
        return v, "Q5"
    if name == "120cell":
        two = Q5(2)
        v = signed_points([zero, zero, two, two], ALL_PERMUTATIONS, zero)
        v += signed_points([one, one, one, sq5], ALL_PERMUTATIONS, zero)
        v += signed_points([iphi * iphi, phi, phi, phi], ALL_PERMUTATIONS, zero)
        v += signed_points([iphi, iphi, iphi, phi * phi], ALL_PERMUTATIONS, zero)
        v += signed_points([zero, iphi * iphi, one, phi * phi], EVEN_PERMUTATIONS, zero)
        v += signed_points([zero, iphi, phi, sq5], EVEN_PERMUTATIONS, zero)
        v += signed_points([iphi, one, phi, two], EVEN_PERMUTATIONS, zero)
        return v, "Q5"
    raise ValueError(f"unknown polytope {name!r}")


# the counts of the literature: vertices, edges, cells, ridges, cell vertices, facets per cell
COUNTS = {
    "5cell": (5, 10, 5, 10, 4, 4),
    "tesseract": (16, 32, 8, 24, 8, 6),
    "16cell": (8, 24, 16, 32, 4, 4),
    "24cell": (24, 96, 24, 96, 6, 8),
    "120cell": (600, 1200, 120, 720, 20, 12),
    "600cell": (120, 720, 600, 1200, 4, 4),
}


def _norm_str(x):
    return x


def dual_normals(name, V):
    """Direction vectors whose maxima over the vertices are the cells (the dual polytope's vertices)."""
    F = Fraction
    if name == "16cell":
        return [[F(s) for s in signs] for signs in itertools.product((-1, 1), repeat=4)]
    if name == "tesseract":
        return [[F(s if i == k else 0) for i in range(4)] for k in range(4) for s in (-1, 1)]
    if name == "24cell":
        n = [[F(s if i == k else 0) for i in range(4)] for k in range(4) for s in (-1, 1)]
        n += [[F(s) for s in signs] for signs in itertools.product((-1, 1), repeat=4)]
        return n
    if name == "120cell":
        # the vertices of the dual 600-cell, in the orientation dual to the coordinates of raw_vertices("120cell"): the
        # same as the 600-cell's with the ODD permutations of (phi, 1, 1/phi, 0) / 2 (checked: every one of the 120
        # directions has a maximum set of 20 vertices)
        zero, one, half = Q5(0), Q5(1), Q5(F(1, 2))
        phi, iphi = Q5(F(1, 2), F(1, 2)), Q5(F(-1, 2), F(1, 2))
        odd = [p for p in ALL_PERMUTATIONS if p not in EVEN_PERMUTATIONS]
        n = signed_points([one, zero, zero, zero], ALL_PERMUTATIONS, zero)
        n += signed_points([half, half, half, half], [(0, 1, 2, 3)], zero)
        n += signed_points([phi * half, half, iphi * half, zero], odd, zero)
        return n
    raise ValueError(name)


class Cell:
    """One cell: exact data and chart. Filled by Container."""


class Container:
    """The exact polytope: vertices, cells in the engine's order, ridges and gluings in charts."""

    def __init__(self, name):
        self.name = name
        V, K = raw_vertices(name)
        # one number type for every container: Q5 with a zero irrational part is a rational with integer arithmetic (faster
        # than Fraction), so the rational containers (5-, 16-, 24-cell, tesseract) and the golden ones run the same code
        V = [[x if isinstance(x, Q5) else Q5(x) for x in v] for v in V]
        self.V, self.K = V, K
        self.n = len(V[0])
        self.zero = Q5(0)
        e = vsub(V[0], V[1])
        # the least squared distance from vertex 0 is the squared edge
        d2 = [vdot(vsub(V[0], w), vsub(V[0], w)) for w in V[1:]]
        self.e2 = d2[0]
        for x in d2:
            if x < self.e2:
                self.e2 = x
        self._adj = None
        self._build_cells()
        self._order_like_engine()
        self._ridges()
        self._charts()

    # -------------------------------------------------------------- combinatorics
    def _edge_graph(self):
        """Adjacency of the edge graph: pairs at exactly the squared edge. Candidates come from floats (1e-9), and each is
        confirmed in exact arithmetic (an exact equality implies a float one within rounding)."""
        nv = len(self.V)
        fV = [[float(x) for x in v] for v in self.V]
        e2f = float(self.e2)
        adj = [[] for _ in range(nv)]
        for i in range(nv):
            for j in range(i + 1, nv):
                d2 = sum((a - b) ** 2 for a, b in zip(fV[i], fV[j]))
                if abs(d2 - e2f) < 1e-9:
                    d = vsub(self.V[i], self.V[j])
                    if vdot(d, d) == self.e2:
                        adj[i].append(j)
                        adj[j].append(i)
        return adj

    def _build_cells(self):
        name, V = self.name, self.V
        nv = len(V)
        adj = self._edge_graph()
        self.adj = adj
        sets = []
        if name == "5cell":
            sets = [[i for i in range(5) if i != j] for j in range(5)]
        elif name == "600cell":
            # cells are the 4-cliques of the edge graph
            nbr = [set(a) for a in adj]
            seen = set()
            for a in range(nv):
                for b in adj[a]:
                    if b < a:
                        continue
                    for c in nbr[a] & nbr[b]:
                        if c < b:
                            continue
                        for d in nbr[a] & nbr[b] & nbr[c]:
                            if d < c:
                                continue
                            seen.add((a, b, c, d))
            sets = [list(t) for t in sorted(seen)]
        else:
            fV = [[float(x) for x in v] for v in V]
            for N in dual_normals(name, V):
                fN = [float(x) for x in N]
                fd = [sum(a * b for a, b in zip(fN, v)) for v in fV]
                top = max(fd)
                # candidates within 1e-9 of the float maximum; the exact maximum and equality decide
                cand = [i for i, x in enumerate(fd) if x > top - 1e-9]
                ex = {i: vdot(N, V[i]) for i in cand}
                m = ex[cand[0]]
                for x in ex.values():
                    if x > m:
                        m = x
                sets.append([i for i in cand if ex[i] == m])
        nV, nE, nC, nR, cv, cf = COUNTS[name]
        if nv != nV or sum(len(a) for a in adj) // 2 != nE or len(sets) != nC or any(len(s) != cv for s in sets):
            raise AssertionError(f"{name}: counts differ from the literature")
        self.raw_cells = [sorted(s) for s in sets]

    def _order_like_engine(self):
        """Orders the cells like packages/boundary/src/polytope.ts (sorted by the outward unit normal with
        compareNormals, normals from the unit-edge float vertices) and computes the engine's origin and frame."""
        fv = float_vertices(self.name, self)
        self.fvertices = fv

        def fnormal(ids):
            # outward unit normal of the cell: orthogonal complement of the edge vectors, pointing away from the centre
            p0 = fv[ids[0]]
            vecs = [[fv[i][k] - p0[k] for k in range(4)] for i in ids[1:]]
            basis = orthonormal(vecs, 1e-9)
            n = orth_complement(basis)
            cen = [sum(fv[i][k] for i in ids) / len(ids) for k in range(4)]
            if sum(n[k] * cen[k] for k in range(4)) < 0:
                n = [-x for x in n]
            return [0.0 if abs(x) < 1e-15 else x for x in n]

        normals = [fnormal(ids) for ids in self.raw_cells]
        order = sorted(range(len(normals)), key=__import__("functools").cmp_to_key(lambda a, b: compare_normals(normals[a], normals[b])))
        self.cells = []
        for idx, ci in enumerate(order):
            c = Cell()
            c.index = idx
            c.ids = self.raw_cells[ci]
            c.fnormal = normals[ci]
            cen = [sum(fv[i][k] for i in c.ids) / len(c.ids) for k in range(4)]
            c.forigin = [0.0 if abs(x) < 1e-15 else x for x in cen]
            std = [[1.0 if i == j else 0.0 for i in range(4)] for j in range(4)]
            rej = [reject(e, [c.fnormal]) for e in std]
            c.fframe = [[0.0 if abs(x) < 1e-15 else x for x in f] for f in orthonormal(rej, 1e-6)[:3]]
            self.cells.append(c)
        self.cell_of_ids = {tuple(c.ids): c.index for c in self.cells}

    def _ridges(self):
        nC = len(self.cells)
        by_vertex = {}
        for c in self.cells:
            for v in c.ids:
                by_vertex.setdefault(v, []).append(c.index)
        shared = {}
        for v, cs in by_vertex.items():
            for a in cs:
                for b in cs:
                    if a < b:
                        shared.setdefault((a, b), []).append(v)
        ridges = {}
        for (a, b), vs in shared.items():
            if len(vs) >= 3:
                ridges[(a, b)] = sorted(vs)
        for c in self.cells:
            c.facets = []
        cnt = {}
        for (a, b), vs in ridges.items():
            poly = self._cyclic(vs)
            for x, y in ((a, b), (b, a)):
                self.cells[x].facets.append({"nbr": y, "ids": poly})
                cnt[x] = cnt.get(x, 0) + 1
        nV, nE, nC_, nR, cv, cf = COUNTS[self.name]
        if len(ridges) != nR or any(cnt.get(c.index, 0) != cf for c in self.cells):
            raise AssertionError(f"{self.name}: ridges differ from the literature")
        for c in self.cells:
            c.facets.sort(key=lambda f: f["nbr"])
        for c in self.cells:
            # matching facet in the neighbour (after every cell's facets are in their final order)
            for f in c.facets:
                f["back"] = next(i for i, g in enumerate(self.cells[f["nbr"]].facets) if g["nbr"] == c.index)
            hinges = set()
            for f in c.facets:
                p = f["ids"]
                for i in range(len(p)):
                    a, b = p[i], p[(i + 1) % len(p)]
                    hinges.add((min(a, b), max(a, b)))
            c.hinges = sorted(hinges)

    def _cyclic(self, vs):
        """The vertices of a ridge (a convex polygon) in cyclic order, by the edge graph."""
        vs = list(vs)
        if len(vs) == 3:
            return vs
        s = set(vs)
        out = [vs[0]]
        prev = None
        while len(out) < len(vs):
            cur = out[-1]
            nxt = [w for w in self.adj[cur] if w in s and w != prev and w not in out]
            if not nxt:
                raise AssertionError("ridge is not a cycle")
            prev = cur
            out.append(nxt[0])
        return out

    # -------------------------------------------------------------- exact geometry
    def _charts(self):
        V, zero = self.V, self.zero
        for c in self.cells:
            P = [V[i] for i in c.ids]
            c.N = [zero] * self.n
            for p in P:
                c.N = vadd(c.N, p)
            c.b = vdot(c.N, P[0])
        # all cells have the same |N|^2
        self.a = vdot(self.cells[0].N, self.cells[0].N)
        for c in self.cells:
            c.O_idx = c.ids[0]
            c.O = V[c.O_idx]
            # three edge vectors from O that make an independent triple
            cand = [i for i in c.ids[1:]]
            best = None
            for tri in itertools.combinations(cand, 3):
                B = [vsub(V[i], c.O) for i in tri]
                G = [[vdot(B[i], B[j]) for j in range(3)] for i in range(3)]
                if sgn(det3(G)) != 0:
                    best = (B, G)
                    break
            c.B, c.G = best
            c.Ginv = inv3(c.G)
            c.chart_v = {i: self.chart(c, V[i]) for i in c.ids}
            cen = [zero] * 3
            for i in c.ids:
                cen = vadd(cen, c.chart_v[i])
            c.centroid = [x / len(c.ids) if not isinstance(x, Q5) else x / Fraction(len(c.ids)) for x in cen]
            for f in c.facets:
                f["chart"] = [c.chart_v[i] for i in f["ids"]]

    def ginv(self, c, fi):
        """The gluing map of facet fi of cell c (computed on first use)."""
        f = c.facets[fi]
        if "ginv" not in f:
            f["ginv"] = self._ginv(c, f)
        return f["ginv"]

    def chart(self, c, x):
        """Chart coordinates of a point x of the hyperplane of c."""
        w = vsub(x, c.O)
        return matvec(c.Ginv, [vdot(b, w) for b in c.B])

    def chart_lin(self, c, w):
        return matvec(c.Ginv, [vdot(b, w) for b in c.B])

    def point(self, c, lam):
        """The point of R^n with chart coordinates lam in c."""
        x = list(c.O)
        for i in range(3):
            x = vadd(x, vscale(lam[i], c.B[i]))
        return x

    def rot(self, N, N2):
        """The rotation (as a function of vectors) taking N to N2 in their plane and fixing the rest."""
        a = self.a
        k = vdot(N, N2) / a
        E = vsub(N2, vscale(k, N))
        one_minus = 1 - k * k
        z = (k - 1) / (one_minus * a) if not isinstance(one_minus, Q5) else (k - 1) / (one_minus * a)

        def r(x):
            nx = vdot(N, x)
            ex = vdot(E, x)
            t1 = vscale((k - 1) / a * nx, N)
            t2 = vscale(z * ex, E)
            t3 = vscale(1 / a if not isinstance(a, Q5) else a.inv(), vsub(vscale(nx, E), vscale(ex, N)))
            return vadd(vadd(vadd(x, t1), t2), t3)

        return r

    def _ginv(self, c, f):
        """The affine map from the chart of the neighbour to the chart of c that develops the neighbour across the
        ridge f into c's hyperplane: lam' = A lam + t."""
        c2 = self.cells[f["nbr"]]
        q = self.V[f["ids"][0]]
        R = self.rot(c2.N, c.N)

        def dev(x):
            return vadd(q, R(vsub(x, q)))

        t = self.chart(c, dev(c2.O))
        cols = [self.chart_lin(c, R(b)) for b in c2.B]
        A = [[cols[j][i] for j in range(3)] for i in range(3)]
        return (A, t)

    # -------------------------------------------------------------- reading a configuration
    def ambient_float(self, cell, u, s):
        c = self.cells[cell]
        x = [s * o for o in c.forigin]
        for k in range(3):
            for i in range(4):
                x[i] += u[k] * c.fframe[k][i]
        return x


_CACHE = {}


def get_container(name):
    if name not in _CACHE:
        _CACHE[name] = Container(name)
    return _CACHE[name]


# ------------------------------------------------------------------ float conventions of the engine


def float_vertices(name, cont=None):
    """The unit-edge float vertices of packages/boundary/src/containers.ts (twin computation in floats)."""
    V, K = raw_vertices(name) if cont is None else (cont.V, cont.K)
    if name == "tesseract":
        return [[float(x) for x in v] for v in V]
    if name == "5cell":
        a = 1 / math.sqrt(5)
        raw = [[1, 1, 1, -a], [1, -1, -1, -a], [-1, 1, -1, -a], [-1, -1, 1, -a], [0, 0, 0, 4 * a]]
        e = math.sqrt(8)
    else:
        raw = [[float(x) for x in v] for v in V]
        e = {"16cell": math.sqrt(2), "24cell": math.sqrt(2), "600cell": 2 / (1 + math.sqrt(5)), "120cell": 3 - math.sqrt(5)}[name]
    return [[x / e for x in v] for v in raw]


def dotf(a, b):
    return sum(x * y for x, y in zip(a, b))


def orthonormal(vs, eps):
    out = []
    for v in vs:
        r = list(v)
        for _ in range(2):
            for b in out:
                c = dotf(r, b)
                r = [x - c * y for x, y in zip(r, b)]
        n = math.sqrt(dotf(r, r))
        if n > eps * max(1.0, math.sqrt(dotf(v, v))):
            out.append([x / n for x in r])
    return out


def reject(v, basis):
    r = list(v)
    for _ in range(2):
        for b in basis:
            c = dotf(r, b)
            r = [x - c * y for x, y in zip(r, b)]
    return r


def orth_complement(basis):
    """The unit vector of R^4 orthogonal to three orthonormal vectors (sign unspecified)."""
    best = None
    for k in range(4):
        e = [1.0 if i == k else 0.0 for i in range(4)]
        r = reject(e, basis)
        n = math.sqrt(dotf(r, r))
        if best is None or n > best[0]:
            best = (n, r)
    return [x / best[0] for x in best[1]]


ORDER_TOL = 1e-9


def compare_normals(a, b):
    """packages/boundary/src/polytope.ts compareNormals."""

    def dominant(n):
        m = max(abs(x) for x in n)
        return next(i for i, x in enumerate(n) if abs(x) >= m - ORDER_TOL)

    ka, kb = dominant(a), dominant(b)
    if ka != kb:
        return ka - kb
    sa = (a[ka] > 0) - (a[ka] < 0)
    sb = (b[kb] > 0) - (b[kb] < 0)
    if sa != sb:
        return sa - sb
    for i in range(len(a)):
        if abs(a[i] - b[i]) > ORDER_TOL:
            return -1 if a[i] < b[i] else 1
    return 0


# ------------------------------------------------------------------ reading a configuration exactly

SNAP = Fraction(1, 10 ** 12)  # a centre outside its closed cell by less than this (times s) is moved onto the cell
WEIGHT_BITS = 60


def _solve(A, y):
    """Least-squares-free solve of A w = y for a square float system by Gauss-Jordan with partial pivoting."""
    n = len(A)
    M = [list(A[i]) + [y[i]] for i in range(n)]
    for c in range(n):
        p = max(range(c, n), key=lambda r: abs(M[r][c]))
        M[c], M[p] = M[p], M[c]
        if abs(M[c][c]) < 1e-12:
            raise ValueError("degenerate")
        for r in range(n):
            if r != c:
                f = M[r][c] / M[c][c]
                M[r] = [a - f * b for a, b in zip(M[r], M[c])]
    return [M[i][n] / M[i][i] for i in range(n)]


def barycentric_weights(cont, cell, x_unit):
    """Dyadic weights w (summing to 1) with sum w_i v_i = x_unit for four independent vertices of the cell (float
    unit-edge vertices of the engine). Returns (vertex ids, weights as Fractions)."""
    c = cont.cells[cell]
    fv = cont.fvertices
    best = None
    for quad in itertools.combinations(c.ids, 4):
        # unknowns w_0..w_3; equations sum w = 1 and the coordinates 0, 1, 2 (the 4th is redundant up to rounding)
        for drop in range(4):
            rows = [k for k in range(4) if k != drop]
            A = [[1.0] * 4] + [[fv[i][k] for i in quad] for k in rows]
            try:
                w = _solve(A, [1.0] + [x_unit[k] for k in rows])
            except ValueError:
                continue
            # residual in the dropped coordinate decides the quality
            res = abs(sum(wi * fv[i][drop] for wi, i in zip(w, quad)) - x_unit[drop])
            if best is None or res < best[0]:
                best = (res, quad, w)
            if res < 1e-13:
                break
        if best and best[0] < 1e-13:
            break
    if best is None:
        raise ValueError("no independent vertices")
    res, quad, w = best
    if res > 1e-9:
        raise ValueError("the point is not on the cell's hyperplane")
    q = Fraction(1, 1 << WEIGHT_BITS)
    ws = [Fraction(round(x * (1 << WEIGHT_BITS)), 1 << WEIGHT_BITS) for x in w]
    big = max(range(4), key=lambda i: abs(ws[i]))
    ws[big] += 1 - sum(ws)
    return list(quad), ws


def facet_planes(cont, c):
    """For each facet of c: (covariant normal pointing out of the cell, in chart coordinates, and a point of the facet)."""
    out = []
    for f in c.facets:
        Q = f["chart"]
        N = cross(vsub(Q[1], Q[0]), vsub(Q[2], Q[0]))
        # orient out of the cell: the centroid is on the negative side
        if sgn(vdot(N, vsub(c.centroid, Q[0]))) > 0:
            N = [-x for x in N]
        out.append((N, Q[0]))
    return out


class Ball:
    """A ball's centre: exact point, and its chart in every closed cell that contains it."""


class Problem:
    """A parsed configuration of balls on a container, all exact."""

    def __init__(self, cont, config, points=None):
        """`points`: exact centres (vectors of R^n in the base polytope's coordinates, one per piece), which replace the reading
        of the floats (for exact constructions and tests; each must lie in the hyperplane of the piece's cell)."""
        self.cont = cont
        s = config["s"]
        self.s = s
        cont_e2 = cont.e2
        # base squared diameter of a ball: e0^2 / s^2
        self.D2 = cont_e2 / (s * s)
        self.balls = []
        zero = cont.zero
        for n, p in enumerate(config["pieces"]):
            cell, u = p["cell"], p["u"]
            b = Ball()
            b.index, b.cell = n, cell
            if points is not None:
                x = list(points[n])
                lam = cont.chart(cont.cells[cell], x)
            else:
                xf = cont.ambient_float(cell, [float(x) for x in u], float(s))
                xu = [x / float(s) for x in xf]
                ids, ws = barycentric_weights(cont, cell, xu)
                x = [zero] * cont.n
                for i, w in zip(ids, ws):
                    x = vadd(x, vscale(w, cont.V[i]))
                lam = cont.chart(cont.cells[cell], x)
                lam, moved = self._snap(cont, cont.cells[cell], lam, s)
                b.snapped = moved
            b.x = cont.point(cont.cells[cell], lam)
            self.balls.append(b)
        for b in self.balls:
            b.cells = {}
            for c in cont.cells:
                if vdot(c.N, b.x) == c.b:
                    b.cells[c.index] = cont.chart(c, b.x)
            if b.cell not in b.cells:
                b.cells[b.cell] = cont.chart(cont.cells[b.cell], b.x)

    @staticmethod
    def _snap(cont, c, lam, s):
        """Moves a centre that is slightly outside its closed cell onto the cell (header); returns (lam, moved)."""
        planes = facet_planes(cont, c)
        tol_sq = SNAP * SNAP * cont.e2  # (1e-12 of an edge)^2, in squared base lengths
        moved = False
        for _ in range(6):
            worst = None
            for N, q in planes:
                ex = vdot(N, vsub(lam, q))
                if sgn(ex) > 0:
                    worst = (N, q, ex) if worst is None or ex > worst[2] else worst
            if worst is None:
                return lam, moved
            N, q, ex = worst
            m = matvec(c.Ginv, N)
            nn = vdot(N, m)
            dist2 = ex * ex / nn
            if dist2 > tol_sq:
                raise ValueError("a centre lies outside its anchor cell")
            lam = vsub(lam, vscale(ex / nn, m))
            moved = True
        raise ValueError("a centre lies outside its anchor cell")


# ------------------------------------------------------------------ chains


class Node:
    __slots__ = ("cell", "A", "t", "window", "gate", "gd2", "entry", "depth", "verts")

    def __init__(self, cell, A, t, window, gate, gd2, entry, depth):
        self.cell, self.A, self.t, self.window, self.gate, self.gd2, self.entry, self.depth = cell, A, t, window, gate, gd2, entry, depth
        self.verts = None

    def dev(self, lam):
        return vadd(matvec(self.A, lam), self.t)


class ChainLimit(Exception):
    pass


def clip_polygon(P, N, u):
    """Sutherland-Hodgman: the part of the planar convex polygon P with N.(x - u) >= 0, exact."""
    d = [vdot(N, vsub(p, u)) for p in P]
    sg = [sgn(x) for x in d]
    if all(x >= 0 for x in sg):
        return P
    if all(x <= 0 for x in sg):
        return [p for p, x in zip(P, sg) if x == 0] if sum(1 for x in sg if x == 0) >= 3 else []
    out = []
    n = len(P)
    for i in range(n):
        j = (i + 1) % n
        if sg[i] >= 0:
            out.append(P[i])
        if (sg[i] > 0 and sg[j] < 0) or (sg[i] < 0 and sg[j] > 0):
            t = d[i] / (d[i] - d[j])
            out.append(vadd(P[i], vscale(t, vsub(P[j], P[i]))))
    # drop repeated neighbours
    ded = []
    for p in out:
        if not ded or ded[-1] != p:
            ded.append(p)
    while len(ded) > 1 and ded[0] == ded[-1]:
        ded.pop()
    return ded


def polygon_area_vec(P):
    s = [P[0][0] * 0] * 3
    for i in range(1, len(P) - 1):
        s = vadd(s, cross(vsub(P[i], P[0]), vsub(P[i + 1], P[0])))
    return s


def seg_dist2(G, u, a, b):
    """Squared distance (metric G) from u to the segment [a, b]."""
    e = vsub(b, a)
    ge = matvec(G, e)
    L = vdot(e, ge)
    w = vsub(u, a)
    if sgn(L) == 0:
        r = vsub(u, a)
        return vdot(r, matvec(G, r))
    t = vdot(w, ge) / L
    if sgn(t) < 0:
        t = t * 0
    elif sgn(t - 1) > 0:
        t = t * 0 + 1
    r = vsub(w, vscale(t, e))
    return vdot(r, matvec(G, r))


def poly_dist2(G, Ginv, u, P):
    """Exact squared distance from u to the planar convex polygon P (chart coordinates, metric G)."""
    Np = polygon_area_vec(P)
    if sgn(Np[0]) == 0 and sgn(Np[1]) == 0 and sgn(Np[2]) == 0:
        return min(seg_dist2(G, u, P[i], P[(i + 1) % len(P)]) for i in range(len(P)))
    m = matvec(Ginv, Np)
    nn = vdot(Np, m)
    h = vdot(Np, vsub(u, P[0]))
    plane2 = h * h / nn
    # the foot of the perpendicular
    f = vsub(u, vscale(h / nn, m))
    inside = True
    for i in range(len(P)):
        a, b = P[i], P[(i + 1) % len(P)]
        c = cross(vsub(b, a), vsub(f, a))
        if sgn(vdot(c, Np)) < 0:
            inside = False
            break
    if inside:
        return plane2
    return min(seg_dist2(G, u, P[i], P[(i + 1) % len(P)]) for i in range(len(P)))


def enumerate_chains(cont, anchor, u, R2, max_nodes=200000):
    """All chains from the point u (chart coordinates in the cell `anchor`) whose gate is within sqrt(R2) (base units
    squared) of u. Returns the nodes, the root first. Exact."""
    G0, Ginv0 = cont.cells[anchor].G, cont.cells[anchor].Ginv
    zero = cont.zero
    I = [[zero + 1 if i == j else zero for j in range(3)] for i in range(3)]
    root = Node(anchor, I, [zero, zero, zero], [], [], zero, -1, 0)
    nodes = [root]
    k = 0
    while k < len(nodes):
        node = nodes[k]
        k += 1
        c = cont.cells[node.cell]
        cen = node.dev(c.centroid)
        for fi, f in enumerate(c.facets):
            if fi == node.entry:
                continue
            Q = [node.dev(p) for p in f["chart"]]
            Np = cross(vsub(Q[1], Q[0]), vsub(Q[2], Q[0]))
            if sgn(vdot(Np, vsub(cen, Q[0]))) > 0:
                Np = [-x for x in Np]
            # u strictly on the inner side of the ridge's plane
            if sgn(vdot(Np, vsub(u, Q[0]))) >= 0:
                continue
            gate = Q
            for N in node.window:
                gate = clip_polygon(gate, N, u)
                if len(gate) < 3:
                    break
            if len(gate) < 3:
                continue
            av = polygon_area_vec(gate)
            if sgn(av[0]) == 0 and sgn(av[1]) == 0 and sgn(av[2]) == 0:
                continue
            gd2 = poly_dist2(G0, Ginv0, u, gate)
            if gd2 > R2:
                continue
            # the cone over the gate
            gc = [zero, zero, zero]
            for p in gate:
                gc = vadd(gc, p)
            gc = [x / len(gate) if not isinstance(x, Q5) else x / Fraction(len(gate)) for x in gc]
            rays = [vsub(p, u) for p in gate]
            ref = vsub(gc, u)
            window = []
            for i in range(len(gate)):
                N = cross(rays[i], rays[(i + 1) % len(gate)])
                sN = sgn(vdot(N, ref))
                if sN == 0:
                    if sgn(N[0]) == 0 and sgn(N[1]) == 0 and sgn(N[2]) == 0:
                        continue
                    # a face through the centroid ray: the gate is degenerate seen from u (cannot happen with area > 0)
                    raise AssertionError("degenerate cone")
                window.append(N if sN > 0 else [-x for x in N])
            if len(nodes) >= max_nodes:
                raise ChainLimit(f"more than {max_nodes} chains")
            c2 = f["nbr"]
            A2, t2 = cont.ginv(c, fi)
            A = matmat(node.A, A2)
            t = vadd(matvec(node.A, t2), node.t)
            nodes.append(Node(c2, A, t, window, gate, gd2, f["back"], node.depth + 1))
    return nodes


def in_window(node, u, y):
    for N in node.window:
        if sgn(vdot(N, vsub(y, u))) < 0:
            return False
    return True


def dist2(G, a, b):
    r = vsub(a, b)
    return vdot(r, matvec(G, r))


def node_vertices(cont, node):
    if node.verts is None:
        c = cont.cells[node.cell]
        node.verts = {i: node.dev(c.chart_v[i]) for i in c.ids}
    return node.verts


def clip_segment(window, u, a, b):
    """The parameter interval [t0, t1] of a + t (b - a), t in [0, 1], inside the closed window; None if empty."""
    t0, t1 = None, None
    one = a[0] * 0 + 1
    zero = a[0] * 0
    t0, t1 = zero, one
    for N in window:
        fa = vdot(N, vsub(a, u))
        fb = vdot(N, vsub(b, u))
        sa, sb = sgn(fa), sgn(fb)
        if sa >= 0 and sb >= 0:
            continue
        if sa < 0 and sb < 0:
            return None
        t = fa / (fa - fb)
        if sa < 0:  # entering
            if t > t0:
                t0 = t
        else:
            if t < t1:
                t1 = t
        if t0 > t1:
            return None
    return t0, t1


def clearance_d2(cont, nodes, u, bound2):
    """The least squared distance (base units, metric of the anchor) from u to the hinges over the chains, considering
    only chains whose gate is within bound2; returns None when no hinge is that close."""
    G0 = cont.cells[nodes[0].cell].G
    best = None
    for node in nodes:
        if sgn(node.gd2 - bound2) > 0:
            continue
        V = node_vertices(cont, node)
        for a, b in cont.cells[node.cell].hinges:
            iv = clip_segment(node.window, u, V[a], V[b])
            if iv is None:
                continue
            t0, t1 = iv
            e = vsub(V[b], V[a])
            ge = matvec(G0, e)
            L = vdot(e, ge)
            if sgn(L) == 0:
                continue
            t = vdot(vsub(u, V[a]), ge) / L
            if t < t0:
                t = t0
            elif t > t1:
                t = t1
            p = vadd(V[a], vscale(t, e))
            d2 = dist2(G0, u, p)
            if best is None or d2 < best:
                best = d2
    return best


class Verdict:
    def __init__(self, strict):
        self.strict = strict
        self.reasons = []
        self.margins = {"clearance_sq": None, "distance_sq": None, "loop_sq": None}
        # per ball, per pair and per loop: the least squared site length found (exact numbers of the field)
        self.clear_d2, self.pair_d2, self.loop_d2 = {}, {}, {}
        self.stats = {}

    @property
    def ok(self):
        return not self.reasons

    def fail(self, msg):
        if len(self.reasons) < 50:
            self.reasons.append(msg)


def _bad(value, limit, strict):
    return sgn(value - limit) <= 0 if strict else sgn(value - limit) < 0


def check_balls(polytope, config, strict=True, reach_factor=None, max_nodes=200000, full=False, points=None):
    """Decides validity of balls (diameter 1) on the boundary of `polytope`. config = {s: Fraction, pieces: [{cell, u}]}
    with exact s; u are floats or Fractions. Returns a Verdict; margins are squared SITE lengths (exact numbers of the field):
    the least clearance^2, the least pair distance^2 and the least loop^2 over the chains examined."""
    cont = get_container(polytope)
    v = Verdict(strict)
    try:
        prob = Problem(cont, config, points)
    except ValueError as e:
        v.fail(str(e))
        return v
    s = prob.s
    D2 = prob.D2
    R2 = D2 if reach_factor is None else D2 * reach_factor
    half2 = D2 / 4
    unit2 = s * s / cont.e2  # site squared length per base squared length
    P = len(prob.balls)
    # runs: one per (ball, closed cell containing it)
    per_ball = []
    for b in prob.balls:
        runs = []
        for ci, lam in b.cells.items():
            try:
                nodes = enumerate_chains(cont, ci, lam, R2, max_nodes)
            except ChainLimit as e:
                v.fail(f"ball {b.index}: {e}")
                return v
            runs.append((ci, lam, nodes))
        per_ball.append(runs)
    v.stats["chains"] = sum(len(r[2]) for runs in per_ball for r in runs)
    v.chains_by_ball = [sum(len(r[2]) for r in runs) for runs in per_ball]
    clear_min = None
    for b, runs in zip(prob.balls, per_ball):
        # clearance: along every chain of every run
        for ci, lam, nodes in runs:
            d = clearance_d2(cont, nodes, lam, half2 if reach_factor is None else R2)
            if d is not None and (clear_min is None or d < clear_min):
                clear_min = d
            if d is not None and (b.index not in v.clear_d2 or d * unit2 < v.clear_d2[b.index]):
                v.clear_d2[b.index] = d * unit2
            if d is not None and _bad(d, half2, strict):
                v.fail(f"ball {b.index} is too close to an edge of the {polytope}: it must be {'more than' if strict else 'at least'} 1/2 from every edge")
    if clear_min is not None:
        v.margins["clearance_sq"] = clear_min * unit2
    if not v.ok and not full:
        return v
    pair_min = None
    loop_min = None
    for i, (b, runs) in enumerate(zip(prob.balls, per_ball)):
        for ci, lam, nodes in runs:
            G0 = cont.cells[ci].G
            by_cell = {}
            for node in nodes:
                by_cell.setdefault(node.cell, []).append(node)
            # other balls (higher index) and this ball (loops)
            for j in range(i, P):
                b2 = prob.balls[j]
                for cj, lam2 in b2.cells.items():
                    for node in by_cell.get(cj, ()):
                        y = node.dev(lam2)
                        if j == i:
                            if node.depth == 0 or y == lam:
                                continue
                        if not in_window(node, lam, y):
                            continue
                        d = dist2(G0, y, lam)
                        tab = v.loop_d2 if j == i else v.pair_d2
                        key = i if j == i else (i, j)
                        if key not in tab or d * unit2 < tab[key]:
                            tab[key] = d * unit2
                        if j == i:
                            if loop_min is None or d < loop_min:
                                loop_min = d
                            if _bad(d, D2, strict):
                                v.fail(f"ball {i} overlaps itself: a geodesic loop leaves and returns to its centre and is shorter than 1")
                        else:
                            if pair_min is None or d < pair_min:
                                pair_min = d
                            if _bad(d, D2, strict):
                                v.fail(f"balls {i} and {j} overlap: a path along the boundary joins their centres and is {'at most' if strict else 'less than'} 1 long")
    v.margins["distance_sq"] = None if pair_min is None else pair_min * unit2
    v.margins["loop_sq"] = None if loop_min is None else loop_min * unit2
    return v


# ------------------------------------------------------------------ the front end used by certify.py

MAX_S = Fraction(10) ** 4
PI_LO = Fraction("3.14159265358979323846264338327950288419716939937510")
_S2 = 10 ** 45
SQRT2_LO = Fraction(math.isqrt(2 * _S2 * _S2), _S2)
SQRT2_HI = SQRT2_LO + Fraction(1, _S2)


class GeneralError(Exception):
    """Malformed input (certify.CertError subclass at the front end)."""


def parse_config(polytope, P, config, to_fraction):
    """The configuration as {s: Fraction, pieces: [{cell, u: [Fraction x 3]}]}, checked for form (ball, cell range, number of pieces)."""
    cont = get_container(polytope)
    if not isinstance(config, dict):
        raise GeneralError("the configuration must be an object with s and pieces")
    s = to_fraction(config.get("s"), "s")
    if not 0 < s <= MAX_S:
        raise GeneralError("s must be greater than 0 and at most 1e4")
    pieces = config.get("pieces")
    if not isinstance(pieces, list):
        raise GeneralError("the configuration needs a list of pieces")
    if len(pieces) != P:
        raise GeneralError(f"the configuration has {len(pieces)} pieces; the problem needs {P}")
    out = []
    C = len(cont.cells)
    for n, p in enumerate(pieces):
        if not isinstance(p, dict):
            raise GeneralError(f"piece {n} must be an object with cell and u")
        extra = set(p) - {"cell", "u", "q"}
        if extra:
            raise GeneralError(f"piece {n} has unknown fields {sorted(extra)}")
        if "q" in p:
            raise GeneralError(f"piece {n} is a ball and must not have a rotation q")
        c = p.get("cell")
        if isinstance(c, bool) or not isinstance(c, int) or not 0 <= c < C:
            raise GeneralError(f"piece {n}: cell must be a whole number from 0 to {C - 1}")
        u = p.get("u")
        if not isinstance(u, list) or len(u) != 3:
            raise GeneralError(f"piece {n}: u must be 3 numbers")
        uu = [to_fraction(x, f"piece {n} u[{k}]") for k, x in enumerate(u)]
        if any(abs(x) > MAX_S for x in uu):
            raise GeneralError(f"piece {n}: u is out of range")
        out.append({"cell": c, "u": uu})
    return {"s": s, "pieces": out}


def summary(v):
    """Readable margins of a verdict (site lengths, rounded down so that they are never overstated)."""
    m = v.margins
    def low(x):
        # a binary64 whose square is not above the exact squared length
        if x is None:
            return None
        lo = x.bounds()[0] if isinstance(x, Q5) else Fraction(x)
        f = math.sqrt(float(lo))
        while f > 0 and Fraction(f) ** 2 > lo:
            f = math.nextafter(f, 0.0)
        return f

    out = {"checker": "python/boundary/general.py", "version": "boundary-general 1", "strict": v.strict, "valid": v.ok}
    out["min_clearance"] = low(m["clearance_sq"])
    out["min_distance"] = low(m["distance_sq"])
    out["min_loop"] = low(m["loop_sq"])
    out["chains"] = v.stats.get("chains")
    return out


def exact_bounds(polytope, P):
    """The proven lower bounds on s for P balls (BOUNDARY.md section 5) as (kind, text, check), where check(x) decides
    that the rational x is not above the bound (with the enclosures of pi, sqrt 2 and sqrt 5 to 45 digits)."""
    F = Fraction
    # volume of the boundary at unit edge: (a + b r) with r = sqrt 2 or sqrt 5
    vol = {
        "5cell": (F(0), F(5, 12), "2"),
        "16cell": (F(0), F(4, 3), "2"),
        "24cell": (F(0), F(8), "2"),
        "600cell": (F(0), F(50), "2"),
        "120cell": (F(450), F(210), "5"),
    }
    # rho^2, the squared distance from a cell's centre to its edges
    rho2 = {"5cell": Q5(F(1, 8)), "16cell": Q5(F(1, 8)), "600cell": Q5(F(1, 8)), "24cell": Q5(F(1, 4)), "120cell": Q5(F(7, 8), F(3, 8))}

    def volume_check(x):
        a, b, r = vol[polytope]
        hi = a + b * (SQRT2_HI if r == "2" else SQRT5_HI)
        return 6 * hi * x ** 3 <= P * PI_LO

    def one_check(x):
        return sgn(4 * rho2[polytope] * (x * x) - 1) <= 0

    return [("volume", "(P pi / (6 V))^(1/3)", volume_check), ("one piece", "1 / (2 rho)", one_check)]


def main(argv):
    import json
    import sys

    if len(argv) != 2:
        print("usage: python3 python/boundary/general.py boundary/<polytope>/ball/<P> <config.json>")
        return 2
    import certify

    with open(argv[1]) as fh:
        config = json.load(fh)
    try:
        print(json.dumps(certify.certify(argv[0], config), indent=2))
    except certify.CertError as e:
        print(f"FAIL {e}")
        return 1
    return 0


if __name__ == "__main__":
    import sys

    sys.exit(main(sys.argv[1:]))
