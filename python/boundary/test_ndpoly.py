"""Tests of the nD certifier (ndpoly.py): balls on the boundary of the d-simplex, d-orthoplex and d-demicube, d = 5..10 (BOUNDARY.md 4D).

  * the containers: facet and ridge counts against the literature, every ridge point on both planes, the vertex sets read by brute force
    (exact rank) against the engine-independent combinatorics (neighbours, hinge pairs) for d = 5, 6 (7 for the hinges of a few cells);
  * the development maps: the reflection with normal n1 + n2 is an exact isometry, an involution, fixes the ridge, maps the neighbour
    onto the far side; it is rational on the simplex, orthoplex and the demicube's right-angled ridges and needs sqrt d on the ridge
    between a demicube facet and a simplex facet (d = 5, 6, 7, 8, 10; for d = 9, sqrt d = 3 is rational);
  * Q(sqrt k): field arithmetic and signs;
  * the lower bounds against the engine's (fixtures/nd_bounds.json), both directions (the value passes, the next float up fails);
  * the reading of a centre: the exact point differs from the engine's floats by under 1e-15 s;
  * the engine's distances (fixtures/nd_crosscheck.json, written by packages/boundary/scripts/gen-nd-crosscheck.ts) against the exact ones;
  * adversarial cases (pairs near 1 apart, balls near 1/2 from a hinge, balls on a ridge): the float validator at tolerance 0 never accepts
    what the exact certifier refuses beyond the reading discrepancy, and the exact strict verdict implies the engine's;
  * the test packings of fixtures/nd_packings.json certified strictly valid with margins of 1e-12 and more, shrunk copies refused, the exact
    constructions (simplex d+1 balls at the cells' centres from sqrt(d(d-1)/2), demicube 2d from sqrt 2, orthoplex 2) at 1 +- 1e-9;
  * fail-closed: a step limit or a chain limit raises, never accepts; malformed input; certify.py's entry point.

Run: python3 python/boundary/test_ndpoly.py
"""

import json
import math
import os
import random
import sys
import time
import unittest
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import certify as C  # noqa: E402
import ndpoly as N  # noqa: E402

FIX = os.path.join(HERE, "fixtures")
tf = C.to_fraction


def load(name):
    with open(os.path.join(FIX, name)) as fh:
        return json.load(fh)


def cfg(s, pieces):
    return {"s": s, "pieces": [{"cell": c, "u": list(u)} for c, u in pieces]}


def rank_exact(rows):
    """The rank of a list of rational vectors (Gauss-Jordan over Fractions)."""
    rows = [list(r) for r in rows]
    rk = 0
    ncol = len(rows[0]) if rows else 0
    for col in range(ncol):
        piv = next((i for i in range(rk, len(rows)) if rows[i][col] != 0), None)
        if piv is None:
            continue
        rows[rk], rows[piv] = rows[piv], rows[rk]
        for i in range(len(rows)):
            if i != rk and rows[i][col] != 0:
                f = rows[i][col] / rows[rk][col]
                rows[i] = [a - f * b for a, b in zip(rows[i], rows[rk])]
        rk += 1
    return rk


def affine_rank(vs):
    if len(vs) < 2:
        return 0
    return rank_exact([[a - b for a, b in zip(v, vs[0])] for v in vs[1:]])


class Containers(unittest.TestCase):
    def test_counts(self):
        for d in range(5, 11):
            self.assertEqual(N.get_cont(f"{d}simplex").C, d + 1)
            self.assertEqual(N.get_cont(f"{d}orthoplex").C, 2 ** d)
            self.assertEqual(N.get_cont(f"{d}demicube").C, 2 * d + 2 ** (d - 1))
        # the 5-demicube has 16 + 10 facets
        self.assertEqual(N.get_cont("5demicube").C, 26)
        self.assertEqual(N.get_cont("10orthoplex").C, 1024)
        self.assertEqual(N.get_cont("10demicube").C, 532)

    def test_normals_have_one_squared_length_and_ridge_points_lie_on_both_planes(self):
        for name in N.NAMES:
            cn = N.get_cont(name)
            cells = range(cn.C) if cn.C <= 70 else [0, 1, 5, 17, 33, 77, cn.C - 1]
            for c in cells:
                cl = cn.cell(c)
                self.assertEqual(N.dot(cl.N, cl.N), cn.kappa, (name, c))
                for k, g in enumerate(cl.nbr):
                    q = cl.rp[k]
                    self.assertEqual(N.sg(N.dot(cl.N, q) - cl.h), 0, (name, c, g))
                    self.assertEqual(N.sg(N.dot(cl.Nn[k], q) - cl.hn[k]), 0, (name, c, g))
                    self.assertIn(c, cn.cell(g).nbr)

    def test_brute_force_neighbours_and_hinges_from_the_vertex_set(self):
        for name in [n for n in N.NAMES if n[0] in "56"]:
            cn = N.get_cont(name)
            # all vertices of the container, from the cells' vertex lists
            verts = {}
            for c in range(cn.C):
                for v in cn.vertices(c):
                    verts[tuple(v)] = None
            V = [list(v) for v in verts]
            # every vertex satisfies every facet inequality; each cell's vertices span n dimensions
            cells = range(cn.C) if cn.C <= 70 else [0, 3, 40, 63]
            for c in cells:
                cl = cn.cell(c)
                for v in V:
                    self.assertLessEqual(N.sg(N.dot(cl.N, v) - cl.h), 0)
                tv = [v for v in V if N.sg(N.dot(cl.N, v) - cl.h) == 0]
                self.assertEqual(affine_rank(tv), cn.n, (name, c))
                brute = []
                for g in range(cn.C):
                    if g == c:
                        continue
                    gl = cn.cell(g)
                    common = [v for v in tv if N.sg(N.dot(gl.N, v) - gl.h) == 0]
                    if len(common) >= cn.n and affine_rank(common) == cn.n - 1:
                        brute.append(g)
                self.assertEqual(sorted(cl.nbr), brute, (name, c))
                pairs = set()
                for k in range(len(cl.nbr)):
                    for l in range(k + 1, len(cl.nbr)):
                        gk, gl_ = cn.cell(cl.nbr[k]), cn.cell(cl.nbr[l])
                        common = [v for v in tv if N.sg(N.dot(gk.N, v) - gk.h) == 0 and N.sg(N.dot(gl_.N, v) - gl_.h) == 0]
                        if len(common) >= cn.n - 1 and affine_rank(common) == cn.n - 2:
                            pairs.add((k, l))
                self.assertEqual(set(cl.hinges), pairs, (name, c))

    def test_hinge_pairs_of_the_7_demicube(self):
        cn = N.get_cont("7demicube")
        verts = {}
        for c in range(cn.C):
            for v in cn.vertices(c):
                verts[tuple(v)] = None
        V = [list(v) for v in verts]
        for c in (0, 1, 70):
            cl = cn.cell(c)
            tv = [v for v in V if N.sg(N.dot(cl.N, v) - cl.h) == 0]
            pairs = set()
            for k in range(len(cl.nbr)):
                for l in range(k + 1, len(cl.nbr)):
                    gk, gl_ = cn.cell(cl.nbr[k]), cn.cell(cl.nbr[l])
                    common = [v for v in tv if N.sg(N.dot(gk.N, v) - gk.h) == 0 and N.sg(N.dot(gl_.N, v) - gl_.h) == 0]
                    if len(common) >= cn.n - 1 and affine_rank(common) == cn.n - 2:
                        pairs.add((k, l))
            self.assertEqual(set(cl.hinges), pairs, c)


class Developments(unittest.TestCase):
    def sigma(self, cn, c, k):
        cl = cn.cell(c)
        m, mm = N.mirror(cn, cl, k)
        q = cl.rp[k]
        return lambda x: [xi - 2 * N.dot(m, [a - b for a, b in zip(x, q)]) / mm * mi for xi, mi in zip(x, m)]

    def test_the_reflection_is_an_exact_isometry_that_fixes_the_ridge_and_lays_the_neighbour_on_the_far_side(self):
        rnd = random.Random(3)
        for name in N.NAMES:
            cn = N.get_cont(name)
            for _ in range(3):
                c = rnd.randrange(min(cn.C, 60))
                cl = cn.cell(c)
                k = rnd.randrange(len(cl.nbr))
                g = cn.cell(cl.nbr[k])
                sg_ = self.sigma(cn, c, k)
                a = [Fr(rnd.randint(-5, 5), rnd.randint(1, 4)) for _ in range(cn.D)]
                b = [Fr(rnd.randint(-5, 5), rnd.randint(1, 4)) for _ in range(cn.D)]
                da = [x - y for x, y in zip(sg_(a), sg_(b))]
                dd = [x - y for x, y in zip(a, b)]
                self.assertEqual(N.dot(da, da), N.dot(dd, dd), name)
                self.assertTrue(all(N.sg(x - y) == 0 for x, y in zip(sg_(sg_(a)), a)))
                for v in cn.vertices(c)[:4]:
                    # ridge vertices are fixed: those of cell c that lie on the plane of g
                    if N.sg(N.dot(g.N, v) - g.h) == 0:
                        self.assertTrue(all(N.sg(x - y) == 0 for x, y in zip(sg_(v), v)))
                # the neighbour's centre lands on the plane of c, beyond the ridge: N_g . sigma(o_g) > h_g
                og = [sum(col) / len(cn.vertices(cl.nbr[k])) for col in zip(*cn.vertices(cl.nbr[k]))]
                img = sg_(og)
                self.assertEqual(N.sg(N.dot(cl.N, img) - cl.h), 0)
                self.assertGreater(N.sg(N.dot(g.N, img) - g.h), 0)

    def test_rationality(self):
        # a rational point of the neighbour's plane maps to a rational point exactly when the reflection is rational on it
        def irrational(name, c, k):
            cn = N.get_cont(name)
            sg_ = self.sigma(cn, c, k)
            og = [sum(col) / len(cn.vertices(cn.cell(c).nbr[k])) for col in zip(*cn.vertices(cn.cell(c).nbr[k]))]
            img = sg_(og)
            return any(hasattr(x, "b") and x.b != 0 for x in img)

        for d in range(5, 11):
            for fam in ("simplex", "orthoplex"):
                cn = N.get_cont(f"{d}{fam}")
                for k in range(len(cn.cell(0).nbr)):
                    self.assertFalse(irrational(f"{d}{fam}", 0, k))
            cn = N.get_cont(f"{d}demicube")
            cl = cn.cell(0)  # the facet x_1 = -1/2
            seen = {"dem": False, "sim": False}
            for k, g in enumerate(cl.nbr):
                kind = "dem" if g < 2 * d else "sim"
                seen[kind] = True
                irr = irrational(f"{d}demicube", 0, k)
                if kind == "dem":
                    self.assertFalse(irr, (d, g))
                else:
                    self.assertEqual(irr, d != 9, (d, g))
                if kind == "sim" and k > 40:
                    break
            self.assertTrue(all(seen.values()))


class Fields(unittest.TestCase):
    def test_qf(self):
        for k in (2, 5, 6, 7, 10):
            QF = N.make_qf(k)
            r = QF(0, 1)
            self.assertEqual(r * r, QF(k, 0))
            x = QF(3, -2)
            self.assertEqual(x * x.inv(), QF(1, 0))
            self.assertEqual(x / x, QF(1, 0))
            self.assertEqual(N.sg(r), 1)
            self.assertEqual(N.sg(QF(1, 0) - r * Fr(1, 100)), 1)
            # 3 - 2 sqrt k has the sign of 9 - 4 k
            self.assertEqual(N.sg(x), 1 if 9 > 4 * k else -1)
            self.assertTrue(QF(3, -2) < QF(3, -1))
            lo, hi = r.bounds(30)
            self.assertLess(lo * lo, k)
            self.assertGreater(hi * hi, k)
        self.assertEqual(N.squarefree(8), (2, 2))
        self.assertEqual(N.squarefree(9), (3, 1))
        self.assertEqual(N.squarefree(10), (1, 10))


class Bounds(unittest.TestCase):
    def test_lower_bounds_against_the_engines(self):
        data = load("nd_bounds.json")["bounds"]
        self.assertEqual(len(data), 432)
        for r in data:
            checks = N.exact_bounds(r["polytope"], r["P"])
            x = Fr(r["value"])
            self.assertTrue(any(ok(x) for (_, _, ok) in checks), r)
            up = Fr(math.nextafter(r["value"], math.inf))
            self.assertFalse(any(ok(up) for (_, _, ok) in checks), r)

    def test_one_ball_values(self):
        for d in range(5, 11):
            # simplex, orthoplex: sqrt(d (d-2)) / 2
            for fam in ("simplex", "orthoplex"):
                chk = dict((k, ok) for k, _, ok in N.exact_bounds(f"{d}{fam}", 1))["one piece"]
                b = Fr(math.sqrt(d * (d - 2)) / 2)
                self.assertTrue(chk(b - Fr(1, 10 ** 12)))
                self.assertFalse(chk(b + Fr(1, 10 ** 12)))

    def test_volumes(self):
        # 5-simplex: 6 facets, regular 4-simplices of edge 1, each sqrt 5 / 96
        cn = N.get_cont("5simplex")
        self.assertAlmostEqual(float(N.boundary_volume_hi(cn)), 6 * math.sqrt(5) / 96, places=12)
        cn = N.get_cont("5demicube")
        self.assertAlmostEqual(float(N.boundary_volume_hi(cn)), 10 / 6 + 16 * math.sqrt(5) / 96, places=12)


class Reading(unittest.TestCase):
    def test_the_exact_point_is_within_1e_15_of_the_engines(self):
        rnd = random.Random(1)
        worst = 0.0
        for name in N.NAMES:
            cn = N.get_cont(name)
            for _ in range(3):
                c = rnd.randrange(min(cn.C, 40))
                s = rnd.uniform(0.7, 8)
                u = [rnd.uniform(-0.3, 0.3) * s / math.sqrt(cn.n) for _ in range(cn.n)]
                worst = max(worst, N.reading_error(name, c, u, s))
        self.assertLess(worst, 1e-15)
        # the strict test keeps far more room than that
        self.assertGreater(float(N.READING_CLEAR), 100 * worst)


class EngineCrosscheck(unittest.TestCase):
    def test_engine_distances_equal_the_exact_ones(self):
        data = load("nd_crosscheck.json")["cases"]
        self.assertGreater(len(data), 100)
        worst = 0.0
        compared = 0
        for c in data:
            m = N.measure(c["polytope"], len(c["pieces"]), {"s": c["s"], "pieces": c["pieces"]}, tf)
            sq = lambda x: None if x is None else math.sqrt(float(x))
            for a, b in zip(c["clearance"], m["clearance_sq"]):
                b = sq(b)
                if a is None or b is None:
                    self.assertTrue(a is None and b is None or max(a or 0, b or 0) > 0.4999999999, (c["polytope"], a, b))
                else:
                    worst = max(worst, abs(a - b))
                    compared += 1
                    self.assertLess(abs(a - b), 1e-12, c["polytope"])
            for pr in c["pairs"]:
                b = sq(m["pairs_sq"].get((pr["i"], pr["j"])))
                a = pr["dist"]
                if a is None or b is None:
                    self.assertTrue(a is None and b is None or max(a or 0, b or 0) > 0.999999999999, (c["polytope"], a, b))
                else:
                    worst = max(worst, abs(a - b))
                    compared += 1
                    self.assertLess(abs(a - b), 1e-12, c["polytope"])
            for a, b in zip(c["loops"], m["loops_sq"]):
                self.assertEqual(a is None, b is None, (c["polytope"], a, b))
        print(f"\nengine against exact: {compared} distances compared, largest difference {worst:.2e}")
        self.assertGreater(compared, 150)
        self.assertLess(worst, 1e-12)


class Adversarial(unittest.TestCase):
    """The float validator at tolerance 0 never accepts what the exact certifier refuses (beyond the reading discrepancy 2e-14), and a
    configuration the exact certifier accepts strictly is accepted by the engine."""

    def test_no_false_accept(self):
        data = load("nd_crosscheck.json")["adversarial"]
        self.assertEqual(len(data), 270)

        def dem(c):
            return int(c["polytope"][: -len("demicube")]) if c["polytope"].endswith("demicube") else 0

        # every case on the simplex and orthoplex; on the demicube (slow in exact arithmetic) d <= 6 and, for the ridge cases, d <= 8.
        # All 270 were checked once by hand: engine accepts 54, refuses 216, exact strict valid 54, in 161 s
        chosen = [c for c in data if dem(c) == 0 or dem(c) <= 6 or (c["kind"] == "ridge" and dem(c) <= 8)]
        self.assertGreater(len(chosen), 150)
        t0 = time.time()
        accepted = refused = strict_ok = 0
        worst_accepted = 0.0
        for c in chosen:
            P = len(c["pieces"])
            conf = {"s": c["s"], "pieces": c["pieces"]}
            m = N.measure(c["polytope"], P, conf, tf)
            margins = []
            for x in m["clearance_sq"]:
                if x is not None:
                    margins.append(math.sqrt(float(x)) - 0.5)
            for x in m["pairs_sq"].values():
                if x is not None:
                    margins.append(math.sqrt(float(x)) - 1)
            for x in m["loops_sq"]:
                if x is not None:
                    margins.append(math.sqrt(float(x)) - 1)
            mm = min(margins) if margins else 1.0
            if c["valid0"]:
                accepted += 1
                worst_accepted = min(worst_accepted, mm)
                self.assertGreaterEqual(mm, -2e-14, (c["polytope"], c["kind"], c["delta"]))
            else:
                refused += 1
            v = N.check(c["polytope"], "ball", P, conf, tf, strict=True)
            if v.ok:
                strict_ok += 1
                self.assertTrue(c["valid0"], (c["polytope"], c["kind"], c["delta"], "exact strict verdict valid, engine refuses"))
        print(f"\nadversarial: {len(chosen)} cases, engine accepts {accepted}, refuses {refused}, exact strict valid {strict_ok}; worst exact margin among accepted {worst_accepted:.2e}; {time.time() - t0:.0f} s")
        self.assertGreater(accepted, 20)
        self.assertGreater(refused, 20)


class Coincident(unittest.TestCase):
    """Two different balls at one point overlap (4D.2); only a ball's own zero-length loop (a centre on a ridge) is skipped."""

    def test_stacked_balls_are_refused(self):
        for poly, s in [("5simplex", 6.0), ("5simplex", 2.7), ("7orthoplex", 6.0), ("10orthoplex", 4.48), ("6demicube", 3.0), ("9demicube", 1.0001)]:
            n = int("".join(ch for ch in poly if ch.isdigit())) - 1
            for P in (2, 24):
                conf = {"s": s, "pieces": [{"cell": 0, "u": [0.0] * n} for _ in range(P)]}
                for strict in (False, True):
                    v = N.check(poly, "ball", P, conf, tf, strict=strict)
                    self.assertFalse(v.ok, (poly, s, P, strict))
                with self.assertRaises(C.CertError, msg=(poly, s, P)):
                    C.certify(f"boundary/{poly}/ball/{P}", conf)

    def test_nearly_coincident_balls_are_refused_with_their_distance(self):
        for poly in ("5simplex", "7orthoplex", "6demicube"):
            n = int("".join(ch for ch in poly if ch.isdigit())) - 1
            for off in (0.0, 1e-15, 1e-12, 2e-10, 1e-9, 1e-6):
                u2 = [0.0] * n
                u2[0] = off
                conf = cfg(6.0, [(0, [0.0] * n), (0, u2)])
                self.assertFalse(N.check(poly, "ball", 2, conf, tf, strict=False).ok, (poly, off))
                m = N.measure(poly, 2, conf, tf)
                self.assertEqual(len(m["pairs_sq"]), 1)
                self.assertLess(list(m["pairs_sq"].values())[0], 1, (poly, off))

    def test_the_same_point_under_two_cell_names_is_refused(self):
        for poly in ("6orthoplex", "5simplex", "6demicube"):
            cn = N.get_cont(poly)
            s = 6.0
            c = 0
            for g in sorted({x for x in cn.cell(c).nbr})[:3]:
                common = [v for v in cn.vertices(c) if v in cn.vertices(g)]
                self.assertGreaterEqual(len(common), cn.D - 2)
                xb = [sum(float(v[j]) for v in common) / len(common) for j in range(cn.D)]
                x = [t_ * s / math.sqrt(2) for t_ in xb]
                pcs = []
                for cc in (c, g):
                    o = cn.centre_float(cc)
                    fr = cn.frame_float(cc)
                    pcs.append((cc, [sum(f[j] * (x[j] - s * o[j]) for j in range(cn.D)) for f in fr]))
                self.assertFalse(N.check(poly, "ball", 2, cfg(s, pcs), tf, strict=False).ok, (poly, g))

    def test_a_ball_on_a_ridge_alone_is_not_its_own_overlap(self):
        # the 5-orthoplex's cells 0 and 1 share a ridge; a ball on it (read from either cell) has the trivial loop of length 0 only
        cn = N.get_cont("5orthoplex")
        s = 6.0
        common = [v for v in cn.vertices(0) if v in cn.vertices(1)]
        xb = [sum(float(v[j]) for v in common) / len(common) for j in range(cn.D)]
        x = [t_ * s / math.sqrt(2) for t_ in xb]
        o, fr = cn.centre_float(0), cn.frame_float(0)
        u = [sum(f[j] * (x[j] - s * o[j]) for j in range(cn.D)) for f in fr]
        self.assertTrue(N.check("5orthoplex", "ball", 1, cfg(s, [(0, u)]), tf, strict=False).ok)


class GapCases(unittest.TestCase):
    """The 22 pairs 0.9 apart across a ridge that the float validator once accepted (4D.6): the exact certifier refuses them all."""

    def test_refused(self):
        data = load("nd_gap_cases.json")["cases"]
        self.assertEqual(len(data), 22)
        for c in data:
            v = N.check(c["polytope"], "ball", 2, {"s": c["s"], "pieces": c["pieces"]}, tf, strict=False)
            self.assertFalse(v.ok, (c["label"], c["polytope"], c["s"]))


class Packings(unittest.TestCase):
    def test_test_packings_are_strictly_valid_with_margins_and_shrunk_copies_fail(self):
        data = load("nd_packings.json")["packings"]
        self.assertEqual(len(data), 72)
        for r in data:
            self.assertRegex(r["id"], N.ID_RE)
            v = C.check(r["id"], r["config"], strict=True)
            self.assertTrue(v.ok, (r["id"], v.reasons[:1]))
            summ = N.summary(v)
            # margins of the stored summary: 1e-12 and more where they are visible
            for key, limit in (("min_clearance", 0.5), ("min_distance", 1.0), ("min_loop", 1.0)):
                if r[key] is not None:
                    self.assertGreaterEqual(r[key] - limit, 1e-12 - 1e-15, (r["id"], key))
            # shrunk by 1e-6 each is refused (the 9- and 10-dimensional ones are left out: each takes seconds to refuse)
            if r["polytope"].lstrip("0123456789") != "demicube" or int(r["polytope"][: -len("demicube")]) <= 8:
                f = 1 - 1e-6
                shrunk = {"s": r["config"]["s"] * f, "pieces": [{"cell": p["cell"], "u": [x * f for x in p["u"]]} for p in r["config"]["pieces"]]}
                self.assertFalse(C.check(r["id"], shrunk, strict=False).ok, r["id"])
            del summ

    def test_exact_constructions(self):
        for d in range(5, 11):
            n = d - 1
            z = [0.0] * n
            s0 = math.sqrt(d * (d - 1) / 2)
            # the simplex: d + 1 balls at the cells' centres, valid from sqrt(d (d-1)/2) (adjacent centres 2 r apart across the ridge)
            pieces = [(c, z) for c in range(d + 1)]
            for f, ok in ((1 + 1e-9, True), (1 - 1e-9, False)):
                v = N.check(f"{d}simplex", "ball", d + 1, cfg(s0 * f, pieces), tf, strict=True)
                self.assertEqual(v.ok, ok, (d, f, v.reasons[:1]))
            # the orthoplex: two adjacent cells
            for f, ok in ((1 + 1e-9, True), (1 - 1e-9, False)):
                v = N.check(f"{d}orthoplex", "ball", 2, cfg(s0 * f, [(0, z), (1, z)]), tf, strict=True)
                self.assertEqual(v.ok, ok, (d, f, v.reasons[:1]))
            # the demicube: 2d balls at the centres of the demicube facets, valid from sqrt 2
            if d <= 7:
                pieces = [(c, z) for c in range(2 * d)]
                for f, ok in ((1 + 1e-9, True), (1 - 1e-9, False)):
                    v = N.check(f"{d}demicube", "ball", 2 * d, cfg(math.sqrt(2) * f, pieces), tf, strict=True)
                    self.assertEqual(v.ok, ok, (d, f, v.reasons[:1]))
            # one ball at a cell's centre: from sqrt(d (d-2)) / 2 on the simplex and orthoplex
            b = math.sqrt(d * (d - 2)) / 2
            for fam in ("simplex", "orthoplex"):
                self.assertTrue(N.check(f"{d}{fam}", "ball", 1, cfg(b * (1 + 1e-9), [(0, z)]), tf, strict=True).ok)
                self.assertFalse(N.check(f"{d}{fam}", "ball", 1, cfg(b * (1 - 1e-9), [(0, z)]), tf, strict=True).ok)

    def test_a_ball_on_a_ridge_belongs_to_both_cells(self):
        # a ball centred on the ridge between cells 0 and g of the 5-orthoplex (nudged 1e-13 into cell 0, or out of it, which the exact
        # reading moves back onto the ridge), and a second ball L away in cell g along the common perpendicular: the straight path across the
        # ridge has length L, so L = 0.9 overlaps and L = 1.05 clears, for both nudges
        name, s = "5orthoplex", 3.0
        cn = N.get_cont(name)
        cl = cn.cell(0)
        k = 0
        g = cl.nbr[k]
        gl = cn.cell(g)
        e0 = math.sqrt(2.0)
        q = [float(x) / e0 * s for x in cl.rp[k]]
        kap = float(cn.kappa)
        Nc = [float(x) for x in cl.N]
        Ng = [float(x) for x in gl.N]
        cg = sum(a * b for a, b in zip(Nc, Ng)) / kap
        # the direction in the plane of g, perpendicular to the ridge, pointing into g (away from the plane of cell 0)
        t = [a - cg * b for a, b in zip(Nc, Ng)]
        tn = math.sqrt(sum(x * x for x in t))
        t = [-x / tn for x in t]
        # the in-cell direction at the ridge pointing into cell 0 (away from the plane of g)
        w = [b - cg * a for a, b in zip(Nc, Ng)]
        wn = math.sqrt(sum(x * x for x in w))
        w = [-x / wn for x in w]

        def local(c, x):
            o = cn.centre_float(c)
            fr = cn.frame_float(c)
            return [sum(fr[j][r] * (x[r] - s * o[r]) for r in range(cn.D)) for j in range(cn.n)]

        for nudge in (1e-13, -1e-13):
            xa = [a + nudge * b for a, b in zip(q, w)]
            for L, ok in ((0.9, False), (1.05, True)):
                xb = [a + L * b for a, b in zip(q, t)]
                conf = cfg(s, [(0, local(0, xa)), (g, local(g, xb))])
                v = N.check(name, "ball", 2, conf, tf, strict=False)
                self.assertEqual(v.ok, ok, (nudge, L, v.reasons[:1]))
                # and the straight path across the ridge is found whichever side of it the exact point falls
                if not ok:
                    self.assertTrue(any("overlap" in r for r in v.reasons))


class FailClosedAndInput(unittest.TestCase):
    def test_limits_raise(self):
        z = [0.0] * 5
        # 6-simplex, one ball at a cell's centre at s = 2.6 (bound 2.449): its cell's seven ridges are within 1/2, so 7 chains and a few steps
        one = cfg(2.6, [(0, z)])
        pieces = [(c, z) for c in range(7)]
        many = cfg(math.sqrt(15) * (1 + 1e-9), pieces)
        saved = (N.MAX_STEPS, N.MAX_NODES)
        try:
            N.MAX_STEPS = 3
            with self.assertRaises(N.NDError):
                N.check("6simplex", "ball", 7, many, tf)
            N.MAX_STEPS = saved[0]
            N.MAX_NODES = 2
            with self.assertRaises(N.NDError):
                N.check("6simplex", "ball", 1, one, tf)
        finally:
            N.MAX_STEPS, N.MAX_NODES = saved
        # the same checks pass with the limits restored
        self.assertTrue(N.check("6simplex", "ball", 7, many, tf).ok)
        v = N.check("6simplex", "ball", 1, one, tf)
        self.assertTrue(v.ok)
        self.assertEqual(v.stats["chains"], 7)

    def test_malformed_and_out_of_scope(self):
        z = [0.0] * 4
        good = cfg(3.5, [(0, z)])
        bad = [
            ({"s": 3.5}, "pieces"),
            ({"s": 0, "pieces": good["pieces"]}, "s"),
            ({"s": 1e5, "pieces": good["pieces"]}, "s"),
            ({"s": 3.5, "pieces": [{"cell": 6, "u": z}]}, "cell"),
            ({"s": 3.5, "pieces": [{"cell": -1, "u": z}]}, "cell"),
            ({"s": 3.5, "pieces": [{"cell": True, "u": z}]}, "cell"),
            ({"s": 3.5, "pieces": [{"cell": 0, "u": [0.0] * 3}]}, "u"),
            ({"s": 3.5, "pieces": [{"cell": 0, "u": z, "q": [1, 0, 0, 0]}]}, "rotation"),
            ({"s": 3.5, "pieces": [{"cell": 0, "u": z, "extra": 1}]}, "unknown"),
            ({"s": 3.5, "pieces": [{"cell": 0, "u": [0.0, 0.0, 0.0, 1e9]}]}, "range"),
            ({"s": 3.5, "pieces": [{"cell": 0, "u": [0.0, 0.0, 0.0, float("inf")]}]}, "finite"),
            ({"s": 3.5, "pieces": [{"cell": 0, "u": [0.0, 0.0, 0.0, "x"]}]}, "rational"),
        ]
        for config, what in bad:
            with self.assertRaises(C.CertError, msg=what):
                C.check("boundary/5simplex/ball/1", config)
        with self.assertRaises(C.CertError):
            C.check("boundary/5simplex/ball/2", good)  # one piece for P = 2
        with self.assertRaises(C.CertError):
            C.check("boundary/5simplex/cube/1", good)  # no cubes
        with self.assertRaises(C.CertError):
            C.check("boundary/5simplex/ball/25", good)  # P above 24
        # a centre outside its cell is refused with a reason, by more than 1e-12 of an edge
        out = C.check("boundary/5simplex/ball/1", cfg(3.5, [(0, [10.0, 0, 0, 0])]), strict=False)
        self.assertFalse(out.ok)

    def test_certify_entry_point(self):
        z = [0.0] * 4
        summ = C.certify("boundary/5simplex/ball/1", cfg(1.94, [(0, z)]))
        self.assertTrue(summ["valid"] and summ["checker"] == "python/boundary/ndpoly.py")
        with self.assertRaises(C.CertError):
            C.certify("boundary/5simplex/ball/1", cfg(1.93, [(0, z)]))
        kinds = [k for k, _, _ in C.lower_bounds("ball", 5, "9demicube")]
        self.assertEqual(kinds, ["volume", "one piece"])
        self.assertIsNotNone(C.ID_RE.match("boundary/10orthoplex/ball/24"))
        self.assertIsNone(C.ID_RE.match("boundary/11simplex/ball/1"))
        self.assertIsNone(C.ID_RE.match("boundary/5simplex/ball/01"))
        self.assertIsNone(C.ID_RE.match("boundary/5simplex/ball/1\n"))


if __name__ == "__main__":
    unittest.main(verbosity=1)
