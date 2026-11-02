"""Exact certifier for balls and axis-aligned cubes on the boundary of the d-cube, d = 5..10 ("4D on 5D" to "9D on 10D").

Problem boundary/<d>cube/<ball|cube>/<P> (5cube .. 10cube): P balls of diameter 1 or unit (d-1)-cubes (edge 1, axis-aligned) on
the boundary of the d-cube [-s/2, s/2]^d; the score is the edge s. It is the tesseract's certifier (certify.py, BOUNDARY.md
sections 2 to 4) with 4 replaced by d, written independently of the TypeScript engine, in exact rational arithmetic only:
a binary64 float is a dyadic rational, every number is read exactly and every test is an integer test after putting all
coordinates over one common denominator. Nothing here uses floating point except the readable numbers of the reports.

Conventions (packages/boundary/src/types.ts, BOUNDARY.md 2 and 4B): cell c = 2 i + (sigma > 0) is the facet x_i = sigma h, h = s/2,
i = 0..d-1; its local coordinates u are the other d - 1 coordinates of R^d in increasing order, so the cell is [-h, h]^n,
n = d - 1; a configuration is {s, pieces: [{cell, u}]}; a piece has no rotation (q is refused).

Validity ("fold" placement, BOUNDARY.md 3): (1) each centre lies in its closed anchor cell; (2) no point of the (d-3)-skeleton lies in
the interior of a developed piece; (3) the development is injective on each piece's interior; (4) interiors of distinct pieces
are disjoint. Strict validity: the closed pieces avoid the skeleton, are embedded and pairwise disjoint.

Developments. Every ridge of the d-cube is a right angle, so the cells a straight path crosses unfold into the cubes of the
lattice 2h Z^n; the hinges are the lattice points with at least two odd multiples of h among their coordinates. A path is
described by the lattice planes y_k = (2m + 1) h it crosses in order (a word); each word has an end cell and an exact
development map x = A y + h b, A a signed permutation, b an integer vector, independent of the start point (as in certify.py).
A straight path shorter than 2 s crosses at most two planes of each axis, all in one direction (three consecutive planes are
4h = 2s apart), so the words with |m_k| <= 2 hold every path a ball can need (diameter 1 < 2s as s >= 1/sqrt 2).

Reach. The words are not listed (109,601 per cell for d = 5, 9,864,101 for d = 6, about 10^16 for d = 10) but found for one
pair of points at a time, pruned by the geodesic distance bound, both bounds exact:
  (1) support bound: a path from u whose last plane on axis k is the m_k-th (m_k != 0) is at least
      sqrt(sum_k delta_k^2) long, delta_k = (2|m_k| - 1) h - sign(m_k) u_k, the distance from u to that plane;
  (2) assignment bound: the image of the target is y = 2 h m + yhat with yhat a signed permutation of the target's local
      coordinates, so |y - u|^2 = sum_r (c_r + yhat_r)^2 >= sum_r (|c|_(r) - |u2|_(r))^2 with c = 2 h m - u (both lists
      sorted; rearrangement inequality).
Only the supports (the vectors m) that pass both are expanded into the orders of their crossings, and each ordered word is
tested exactly: its image lies in its closed lattice cube and the segment u -> y crosses its planes at parameters
0 <= t_1 <= t_2 <= ... <= 1. Both bounds hold for every order, so the search is complete for paths up to the reach. `prune=False`
turns both bounds off (a full enumeration, for tests on small cases).
"""

import json
import math
import re
import sys
from fractions import Fraction

NCUBES = {"5cube": 5, "6cube": 6, "7cube": 7, "8cube": 8, "9cube": 9, "10cube": 10}
PIECES = ("ball", "cube")
MAX_P = 9999
MAX_S = Fraction(10) ** 4
ID_RE = re.compile(r"^boundary/(5cube|6cube|7cube|8cube|9cube|10cube)/(ball|cube)/([1-9][0-9]{0,3})\Z")
# pi to 50 decimals, truncated: a rational lower bound
PI_LO = Fraction("3.14159265358979323846264338327950288419716939937510")
# the search budget of one check: steps of the support and order searches
MAX_WORK = 60_000_000
VERSION = "boundary-ncube 1"
CHECKER = "python/boundary/ncube.py"


class NCubeError(Exception):
    """Malformed input, a configuration outside this certifier's scope, or a search over its work limit."""


# ------------------------------------------------------------------ numbers


def sqrt_down(x):
    """The largest binary64 value whose square is <= the rational x >= 0."""
    if x <= 0:
        return 0.0
    f = math.sqrt(float(x))
    while f > 0 and Fraction(f) ** 2 > x:
        f = math.nextafter(f, 0.0)
    while Fraction(math.nextafter(f, math.inf)) ** 2 <= x:
        f = math.nextafter(f, math.inf)
    return f


def float_down(x):
    f = float(x)
    while Fraction(f) > x:
        f = math.nextafter(f, -math.inf)
    while Fraction(math.nextafter(f, math.inf)) <= x:
        f = math.nextafter(f, math.inf)
    return f


def fmt(x):
    return "%.10g" % float(x)


# ------------------------------------------------------------------ the lattice


def axes(d, i):
    return [j for j in range(d) if j != i]


def embed(d, c, u, H):
    """R^d point of local coordinates u of cell c at half-edge H (any exact numbers)."""
    i, sg = c >> 1, (1 if c & 1 else -1)
    x = [0] * d
    x[i] = sg * H
    for k, j in enumerate(axes(d, i)):
        x[j] = u[k]
    return x


def cells_containing(d, x, H):
    """All cells whose closed facet contains the boundary point x (exact)."""
    out = []
    for i in range(d):
        if x[i] == H:
            out.append(2 * i + 1)
        elif x[i] == -H:
            out.append(2 * i)
    return out


def start_map(d, c0):
    """(cell, rows, signs, b): lattice axis k of cell c0 is its k-th local axis; x = A y + h b with A's column k = sign_k e_row_k."""
    i, sg = c0 >> 1, (1 if c0 & 1 else -1)
    rows = list(axes(d, i))
    b = [0] * d
    b[i] = sg
    return (c0, tuple(rows), tuple([1] * (d - 1)), tuple(b))


def step_map(d, w, k, tau):
    """The map after crossing the plane of lattice axis k in direction tau: the unfolding about the ridge, a rotation R of the plane
    (e_i, e_j) with R(sg e_i) = tj e_j and R(tj e_j) = -sg e_i, applied to A and (as a point rotation about the ridge) to b."""
    c, rows, sgn, b = w
    i, sg = c >> 1, (1 if c & 1 else -1)
    j = rows[k]
    tj = tau * sgn[k]
    nrows = list(rows)
    nsgn = list(sgn)
    for r in range(d - 1):
        if rows[r] == j:
            nrows[r] = i
            nsgn[r] = -sg * tj * sgn[r]
    q = [0] * d
    q[i], q[j] = sg, tj
    wv = [b[t] - q[t] for t in range(d)]
    r2 = list(wv)
    r2[i] = -sg * tj * wv[j]
    r2[j] = sg * tj * wv[i]
    nb = tuple(r2[t] + q[t] for t in range(d))
    return (2 * j + (1 if tj > 0 else 0), tuple(nrows), tuple(nsgn), nb)


def clearance_sq(U, H):
    """Squared geodesic distance from the local point U (scaled integers) to the hinges: the two smallest squared slacks
    (H - |U_k|), BOUNDARY.md 4.2: the nearest hinge of the development is an edge-of-cell hinge of the anchor."""
    a = sorted(H - abs(x) for x in U)
    return a[0] * a[0] + a[1] * a[1]


class Work:
    def __init__(self, limit=MAX_WORK):
        self.left = limit
        self.limit = limit

    def spend(self, k=1):
        self.left -= k
        if self.left < 0:
            raise NCubeError(f"the search needed more than {self.limit} steps (too many balls crowd the hinges)")


def realised_sq(d, H, U1, X2, w, m, steps):
    """Squared length (scaled integer) of the straight developed path from U1 to the point X2 along word w with crossings
    `steps` ((k, level), ...), or None when w does not realise it: the image lies in the closed lattice cube m and the segment
    crosses the planes at parameters 0 <= t_1 <= ... <= 1 (a plane the segment is parallel to is not counted)."""
    c, rows, sgn, b = w
    n = d - 1
    Y = []
    for r in range(n):
        y = sgn[r] * (X2[rows[r]] - H * b[rows[r]])
        cen = 2 * H * m[r]
        if y - cen > H or cen - y > H:
            return None
        Y.append(y)
    pn, pd = 0, 1
    for (k, level) in steps:
        den = Y[k] - U1[k]
        if den == 0:
            return None
        num = H * level - U1[k]
        if den < 0:
            num, den = -num, -den
        if num < 0 or num > den:
            return None
        if num * pd < pn * den:
            return None
        pn, pd = num, den
    return sum((Y[r] - U1[r]) ** 2 for r in range(n))


def path_lengths_sq(d, H, c1, U1, X2, reach_sq, exclude_empty=False, prune=True, work=None):
    """Squared lengths (scaled integers) of the realised straight paths from (c1, U1) to the point X2, over all cells that
    contain X2, whose length is at most the reach (squared, scaled); with prune=False every word with |m_k| <= 2 is tried."""
    n = d - 1
    ends = set(cells_containing(d, X2, H))
    if not ends:
        raise NCubeError("a target is not on the boundary")
    work = work or Work()
    e = sorted(abs(x) for x in X2)
    e.remove(H)  # the normal coordinate of the target's facet (its absolute value is H)
    out = []
    m = [0] * n
    start = start_map(d, c1)

    def orders(w, rem, taken, steps, total):
        if len(steps) == total:
            if w[0] in ends:
                L = realised_sq(d, H, U1, X2, w, m, steps)
                if L is not None:
                    out.append(L)
            return
        for k in range(n):
            if rem[k] == 0:
                continue
            work.spend()
            tau = 1 if m[k] > 0 else -1
            rem[k] -= 1
            taken[k] += 1
            orders(step_map(d, w, k, tau), rem, taken, steps + [(k, tau * (2 * (taken[k]) - 1))], total)
            taken[k] -= 1
            rem[k] += 1

    def support(r, cost, total):
        work.spend()
        if r == n:
            if total == 0 and exclude_empty:
                return
            if prune:
                c = sorted(abs(2 * H * m[q] - U1[q]) for q in range(n))
                if sum((c[q] - e[q]) ** 2 for q in range(n)) > reach_sq:
                    return
            orders(start, [abs(x) for x in m], [0] * n, [], total)
            return
        m[r] = 0
        support(r + 1, cost, total)
        for mm in (1, 2):
            for sign in (1, -1):
                dist = (2 * mm - 1) * H - sign * U1[r]
                cp = cost + dist * dist
                if prune and cp > reach_sq:
                    continue
                m[r] = sign * mm
                support(r + 1, cp, total + mm)
        m[r] = 0

    support(0, 0, 0)
    return out


# ------------------------------------------------------------------ configurations


class Cfg:
    """A parsed problem and configuration, all numbers exact."""

    def __init__(self, polytope, piece, P, s, cells, u):
        self.polytope, self.piece, self.P = polytope, piece, P
        self.d = NCUBES[polytope]
        self.n = self.d - 1
        self.s = s
        h = s / 2
        D = math.lcm(2, h.denominator)
        for uu in u:
            for x in uu:
                D = math.lcm(D, x.denominator)
        self.D = D
        self.H = h.numerator * (D // h.denominator)
        self.cells = cells
        self.U = [[x.numerator * (D // x.denominator) for x in uu] for uu in u]
        self.X = [embed(self.d, c, U, self.H) for c, U in zip(cells, self.U)]


def parse_config(polytope, piece, P, config, to_fraction):
    """The exact configuration, checked for form."""
    if polytope not in NCUBES:
        raise NCubeError(f"unknown polytope {polytope!r}")
    if piece not in PIECES:
        raise NCubeError(f"unknown piece {piece!r}")
    if isinstance(P, bool) or not isinstance(P, int) or not 1 <= P <= MAX_P:
        raise NCubeError(f"P must be a whole number from 1 to {MAX_P}")
    d = NCUBES[polytope]
    if not isinstance(config, dict):
        raise NCubeError("the configuration must be an object with s and pieces")
    s = to_fraction(config.get("s"), "s")
    if not 0 < s <= MAX_S:
        raise NCubeError("s must be greater than 0 and at most 1e4")
    pieces = config.get("pieces")
    if not isinstance(pieces, list):
        raise NCubeError("the configuration needs a list of pieces")
    if len(pieces) != P:
        raise NCubeError(f"the configuration has {len(pieces)} pieces; the problem needs {P}")
    cells, us = [], []
    for n_, p in enumerate(pieces):
        if not isinstance(p, dict):
            raise NCubeError(f"piece {n_} must be an object with cell and u")
        extra = set(p) - {"cell", "u"}
        if extra:
            raise NCubeError(f"piece {n_} has unknown fields {sorted(extra)} (a piece on the {polytope} has no rotation q)")
        c = p.get("cell")
        if isinstance(c, bool) or not isinstance(c, int) or not 0 <= c < 2 * d:
            raise NCubeError(f"piece {n_}: cell must be a whole number from 0 to {2 * d - 1}")
        u = p.get("u")
        if not isinstance(u, list) or len(u) != d - 1:
            raise NCubeError(f"piece {n_}: u must be {d - 1} numbers")
        uu = [to_fraction(x, f"piece {n_} u[{k}]") for k, x in enumerate(u)]
        if any(abs(x) > MAX_S for x in uu):
            raise NCubeError(f"piece {n_}: u is out of range")
        cells.append(c)
        us.append(uu)
    return Cfg(polytope, piece, P, s, cells, us)


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


def _bad(value, limit, strict):
    """value < limit (valid needs >=), or value <= limit when strict."""
    return value <= limit if strict else value < limit


def ball_volume_ok(d, P, s):
    """P balls of diameter 1 fit in the boundary's volume 2d s^n? Decided with a lower bound of pi (a failure is certain)."""
    n = d - 1
    num, den, k = {4: (1, 2, 2), 5: (8, 15, 2), 6: (1, 6, 3), 7: (16, 105, 3), 8: (1, 24, 4), 9: (32, 945, 4)}[n]
    vol = Fraction(num, den) * PI_LO ** k / 2 ** n
    return P * vol <= 2 * d * s ** n


def check(polytope, piece, P, config, to_fraction, strict=True, prune=True, work=None):
    cfg = parse_config(polytope, piece, P, config, to_fraction)
    v = Verdict(cfg, strict)
    H = cfg.H
    for n_, U in enumerate(cfg.U):
        if any(abs(x) > H for x in U):
            v.fail(f"piece {n_}: its centre lies outside its anchor cell {cfg.cells[n_]}")
    if not v.ok:
        return v
    if piece == "ball":
        _check_balls(cfg, v, prune, work or Work())
    else:
        _check_cubes(cfg, v)
    return v


def _near_pairs(cfg):
    """Pairs (a, b), a < b, whose centres are at most 1 apart in R^d (exact); the geodesic is never shorter than the chord."""
    D = cfg.D
    r2 = D * D
    # a grid of side 1 on the first four coordinates
    grid = {}
    for n_, X in enumerate(cfg.X):
        grid.setdefault(tuple(x // D for x in X[:4]), []).append(n_)
    out = set()
    import itertools

    for key, members in grid.items():
        for off in itertools.product((-1, 0, 1), repeat=4):
            other = grid.get(tuple(key[t] + off[t] for t in range(4)))
            if not other:
                continue
            for a in members:
                for b in other:
                    if a < b and sum((p - q) ** 2 for p, q in zip(cfg.X[a], cfg.X[b])) <= r2:
                        out.add((a, b))
    return sorted(out)


def _check_balls(cfg, v, prune, work):
    d, H, D, strict = cfg.d, cfg.H, cfg.D, v.strict
    if not ball_volume_ok(d, cfg.P, cfg.s):
        v.fail(f"{cfg.P} balls have more volume than the boundary of the {cfg.polytope} at s = {fmt(cfg.s)}")
        return
    one = D * D
    # the search reaches a little past 1 (0.1 percent), so that margins just above 1 are reported; the verdicts need exactly 1
    reach = (D * 1025) ** 2 // 1024 ** 2
    cl = [clearance_sq(U, H) for U in cfg.U]
    for n_, c in enumerate(cl):
        # clearance >= 1/2: 4 c >= D^2
        if _bad(4 * c, one, strict):
            v.fail(f"ball {n_} is {fmt(sqrt_down(Fraction(c, one)))} from a (d-3)-face of the {cfg.polytope}; it must be {'more than' if strict else 'at least'} 1/2")
    v.margins["clearance_sq"] = Fraction(min(cl), one)
    if not v.ok:
        return
    best_pair = None
    pairs = 0
    for a, b in _near_pairs(cfg):
        pairs += 1
        L = path_lengths_sq(d, H, cfg.cells[a], cfg.U[a], cfg.X[b], reach, prune=prune, work=work)
        if not L:
            continue
        m = min(L)
        best_pair = m if best_pair is None or m < best_pair else best_pair
        if _bad(m, one, strict):
            v.fail(f"balls {a} and {b} overlap: a path of length {fmt(sqrt_down(Fraction(m, one)))} joins their centres; they must be {'more than' if strict else 'at least'} 1 apart")
    best_loop = None
    for a in range(cfg.P):
        L = path_lengths_sq(d, H, cfg.cells[a], cfg.U[a], cfg.X[a], reach, exclude_empty=True, prune=prune, work=work)
        if not L:
            continue
        m = min(L)
        best_loop = m if best_loop is None or m < best_loop else best_loop
        if _bad(m, one, strict):
            v.fail(f"ball {a} overlaps itself: a geodesic loop of length {fmt(sqrt_down(Fraction(m, one)))} leaves and returns to its centre")
    # None: no realised path within the reach (every distance is more than 1.001)
    v.margins["distance_sq"] = None if best_pair is None else Fraction(best_pair, one)
    v.margins["loop_sq"] = None if best_loop is None else Fraction(best_loop, one)
    v.stats["pairs_within_1_in_Rd"] = pairs


def cube_parts(cfg, n_, strict):
    """The parts of axis-aligned cube n as exact boxes (cell, lo, hi) in cell frames, scaled integers. Returns
    (parts, edge_slack, crossing): edge_slack is the second smallest axis slack H - |U_l| - D/2 (valid needs >= 0, strict > 0),
    crossing lists the axes whose faces the cube crosses (slack < 0; strict: <= 0, a face touched from inside counts)."""
    d, H, D = cfg.d, cfg.H, cfg.D
    nn = cfg.n
    half = D // 2
    U = cfg.U[n_]
    c = cfg.cells[n_]
    slack = [H - abs(x) - half for x in U]
    crossing = [k for k in range(nn) if (slack[k] <= 0 if strict else slack[k] < 0)]
    order = sorted(slack)
    lo = [max(U[k] - half, -H) for k in range(nn)]
    hi = [min(U[k] + half, H) for k in range(nn)]
    parts = [(c, tuple(lo), tuple(hi))]
    if len(crossing) == 1:
        k = crossing[0]
        tau = 1 if U[k] > 0 else -1
        assert abs(U[k]) + half <= 3 * H  # s >= 1: the slice lies in the neighbouring lattice cube
        blo, bhi = list(lo), list(hi)
        if tau > 0:
            blo[k], bhi[k] = H, U[k] + half
        else:
            blo[k], bhi[k] = U[k] - half, -H
        cell, rows, sgn, b = step_map(d, start_map(d, c), k, tau)
        ax2 = axes(d, cell >> 1)
        corners = []
        for y in (blo, bhi):
            x = [H * b[t] for t in range(d)]
            for r in range(nn):
                x[rows[r]] += sgn[r] * y[r]
            corners.append([x[a] for a in ax2])
        parts.append((cell, tuple(min(corners[0][t], corners[1][t]) for t in range(nn)), tuple(max(corners[0][t], corners[1][t]) for t in range(nn))))
    return parts, order[1], crossing


def _check_cubes(cfg, v):
    D, strict, d = cfg.D, v.strict, cfg.d
    if cfg.P > 2 * d * cfg.s ** cfg.n:
        v.fail(f"{cfg.P} cubes have more volume than the boundary of the {cfg.polytope} at s = {fmt(cfg.s)}")
        return
    by_cell = {}
    worst_edge = None
    for n_ in range(cfg.P):
        parts, edge, crossing = cube_parts(cfg, n_, strict)
        worst_edge = edge if worst_edge is None or edge < worst_edge else worst_edge
        if len(crossing) >= 2:
            v.fail(f"cube {n_} crosses the faces of its cell in {len(crossing)} axes, so a (d-3)-face of the {cfg.polytope} passes through {'it' if strict else 'its interior'}")
            continue
        for p in parts:
            by_cell.setdefault(p[0], []).append((n_,) + p)
    v.margins["edge_slack"] = Fraction(worst_edge, D)
    if not v.ok:
        return
    worst_sep = None
    nn = cfg.n
    for cell, lst in by_cell.items():
        lst.sort(key=lambda t: t[2][0])
        for x in range(len(lst)):
            a = lst[x]
            for y in range(x + 1, len(lst)):
                b = lst[y]
                g = b[2][0] - a[3][0]
                if g > 0 and worst_sep is not None and g >= worst_sep:
                    break
                sep = max(max(b[2][t] - a[3][t], a[2][t] - b[3][t]) for t in range(nn))
                if worst_sep is None or sep < worst_sep:
                    worst_sep = sep
                if _bad(sep, 0, strict):
                    who = f"cube {a[0]} overlaps itself" if a[0] == b[0] else f"cubes {min(a[0], b[0])} and {max(a[0], b[0])} overlap"
                    v.fail(f"{who} in cell {cell} ({'touching' if sep == 0 else 'depth ' + fmt(Fraction(-sep, D))})")
    v.margins["part_separation"] = None if worst_sep is None else Fraction(worst_sep, D)
    v.stats["parts"] = sum(len(x) for x in by_cell.values())


def summary(v):
    """Readable margins of a verdict (site lengths, rounded down so that they are never overstated)."""
    out = {"checker": CHECKER, "version": VERSION, "strict": v.strict, "valid": v.ok}
    m = v.margins
    if v.cfg.piece == "ball":
        out["min_clearance"] = sqrt_down(m["clearance_sq"]) if "clearance_sq" in m else None
        out["min_distance"] = sqrt_down(m["distance_sq"]) if m.get("distance_sq") is not None else None
        out["min_loop"] = sqrt_down(m["loop_sq"]) if m.get("loop_sq") is not None else None
        out["pairs_within_1_in_Rd"] = v.stats.get("pairs_within_1_in_Rd")
    else:
        out["min_edge_slack"] = float_down(m["edge_slack"]) if "edge_slack" in m else None
        out["min_part_separation"] = float_down(m["part_separation"]) if m.get("part_separation") is not None else None
        out["parts"] = v.stats.get("parts")
    return out


def exact_bounds(polytope, piece, P):
    """The proven lower bounds on s (BOUNDARY.md 4B) as (kind, text, check), check(x) deciding exactly that the rational x is not above the bound."""
    d = NCUBES[polytope]
    n = d - 1
    one = ("one piece", "1/sqrt 2", lambda x: 2 * x * x <= 1) if piece == "ball" else ("one piece", "1", lambda x: x <= 1)
    if piece == "cube":
        vol = ("volume", "(P/2d)^(1/(d-1))", lambda x: 2 * d * x ** n <= P)
    else:
        num, den, k = {4: (1, 2, 2), 5: (8, 15, 2), 6: (1, 6, 3), 7: (16, 105, 3), 8: (1, 24, 4), 9: (32, 945, 4)}[n]
        vol = ("volume", "(P V_(d-1) / (2^d d))^(1/(d-1))", lambda x: 2 ** (n + 1) * d * x ** n <= P * Fraction(num, den) * PI_LO ** k)
    return [vol, one]


def main(argv):
    if len(argv) != 2:
        print("usage: python3 python/boundary/ncube.py boundary/<d>cube/<ball|cube>/<P> <config.json>")
        return 2
    m = ID_RE.match(argv[0])
    if not m:
        print(f"FAIL {argv[0]!r} is not a d-cube problem id")
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
