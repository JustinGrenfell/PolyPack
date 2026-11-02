"""Tests of the d-cube certifier (ncube.py): balls and axis-aligned cubes on the boundary of the d-cube, d = 5..10.

  * the lattice words: the counts of the research note (109,601 per cell for d = 5 with two crossings per axis, 633 and 6,331
    with one), the closed form of the counts, every development map an unfolding (the cube of the word lands on its end cell),
    and the loop lemma (a word that ends in its start cell crosses two axes twice);
  * the reach search against a full enumeration of the words (exact, random pairs, d = 5 and 6), and the clearance closed form
    against the hinge lattice points;
  * the test packings of fixtures/ncube_packings.json certified strictly valid, their shrunk copies refused, the exact tilings
    P = 2d k^(d-1) cubes at s = k (touching: valid, not strictly), ball grids, and the lower bounds against the engine's;
  * adversarial cases: a ball over a (d-3)-face, balls 1e-9 too close across a ridge, balls overlapping around a hinge through
    a word that crosses two ridges, a cube over a hinge, cubes overlapping across a ridge, touching variants of each;
  * the engine's verdicts (fixtures/ncube_crosscheck.json, written by packages/boundary/scripts/gen-ncube-crosscheck.ts);
  * malformed and out-of-scope input, and the work limit.

Run: python3 python/boundary/test_ncube.py
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
import ncube as N  # noqa: E402

FIX = os.path.join(HERE, "fixtures")


def load(name):
    with open(os.path.join(FIX, name)) as fh:
        return json.load(fh)


def words(d, c0, mm):
    """Every word from c0 with at most mm crossings per lattice axis, all in one direction per axis: (map, m, steps)."""
    n = d - 1
    out = []
    stack = [(N.start_map(d, c0), (0,) * n, ())]
    while stack:
        w, m, steps = stack.pop()
        out.append((w, m, steps))
        for k in range(n):
            for tau in (-1, 1):
                if m[k] != 0 and (m[k] > 0) != (tau > 0):
                    continue
                if abs(m[k] + tau) > mm:
                    continue
                m2 = list(m)
                m2[k] += tau
                stack.append((N.step_map(d, w, k, tau), tuple(m2), steps + ((k, 2 * m[k] + tau),)))
    return out


def count_formula(n, mm):
    """The number of words: sequences of distinct crossings, each axis at most mm times in one direction."""
    if mm == 1:
        return sum(math.comb(n, k) * 2 ** k * math.factorial(k) for k in range(n + 1))
    # mm = 2: the exponential generating function (1 + x)^(2n)
    return sum(math.factorial(2 * n) // math.factorial(2 * n - k) for k in range(2 * n + 1))


def num(x):
    """A JSON-able number: Fractions go in as exact strings ("3/2"), as the certifier reads them."""
    return str(x) if isinstance(x, Fr) else x


def cfg(s, pieces):
    return {"s": num(s), "pieces": [{"cell": c, "u": [num(x) for x in u]} for c, u in pieces]}


def pid(d, piece, P):
    return f"boundary/{d}cube/{piece}/{P}"


def verdict(d, piece, config, strict=False, **kw):
    return N.check(f"{d}cube", piece, len(config["pieces"]), config, C.to_fraction, strict=strict, **kw)


class Words(unittest.TestCase):
    def test_counts(self):
        # the research note's counts per cell: 1,957 (3D), 109,601 (4D on 5D), 9,864,101 (5D on 6D); 79, 633, 6,331 with one crossing
        self.assertEqual([count_formula(n, 2) for n in (3, 4, 5)], [1957, 109601, 9864101])
        self.assertEqual([count_formula(n, 1) for n in (3, 4, 5)], [79, 633, 6331])
        self.assertEqual(len(words(5, 0, 2)), 109601)
        self.assertEqual(len(words(5, 3, 1)), 633)
        self.assertEqual(len(words(6, 7, 1)), 6331)
        self.assertEqual(len(words(7, 1, 1)), count_formula(6, 1))
        # no cell is special
        self.assertEqual(len(words(5, 9, 2)), 109601)

    def test_maps_are_unfoldings(self):
        # the cube of a word develops onto its end cell (its centre y = 2 h m lands on the centre of the end facet), the columns of
        # A are a signed permutation tangent to the end cell, and a word and its one-step extension agree on the ridge they share
        for d, c0, mm in ((5, 0, 2), (5, 7, 1), (6, 4, 1), (7, 13, 1)):
            n = d - 1
            ws = {steps: (w, m) for (w, m, steps) in words(d, c0, mm)}
            for steps, (w, m) in ws.items():
                cell, rows, sgn, b = w
                i, sg = cell >> 1, (1 if cell & 1 else -1)
                self.assertEqual(len(set(rows)), n)
                self.assertNotIn(i, rows)
                self.assertTrue(all(x in (-1, 1) for x in sgn))
                centre = [b[t] for t in range(d)]
                for r in range(n):
                    centre[rows[r]] += sgn[r] * 2 * m[r]
                self.assertEqual(centre, [sg if t == i else 0 for t in range(d)], (d, steps))
                if not steps:
                    continue
                parent, _ = ws[steps[:-1]]
                k, level = steps[-1]
                # the columns other than k are those of the parent
                for r in range(n):
                    if r != k:
                        self.assertEqual((w[1][r], w[2][r]), (parent[1][r], parent[2][r]))
                # a point of the crossed plane maps alike by both
                face = [0] * n
                face[k] = level
                for r in range(n):
                    if r != k:
                        face[r] = (r % 3) - 1
                def img(map_):
                    x = [map_[3][t] for t in range(d)]
                    for r in range(n):
                        x[map_[1][r]] += map_[2][r] * face[r]
                    return x
                self.assertEqual(img(w), img(parent), (d, steps))

    def test_holonomy_around_a_hinge(self):
        # two straight ways round the hinge u0 = u1 = h of a cell reach the same cell with maps that differ by a quarter turn
        # about the hinge (three cells round a (d-3)-face: 270 degrees)
        for d in (5, 6, 9):
            n = d - 1
            start = N.start_map(d, 0)
            a = N.step_map(d, start, 1, 1)
            b = N.step_map(d, N.step_map(d, start, 0, 1), 1, 1)
            self.assertEqual(a[0], b[0])
            self.assertNotEqual((a[1], a[2]), (b[1], b[2]))
            # the hinge: y0 = y1 = h (levels at h = 1), the other coordinates free: its points map alike
            for t in (-1, 0, 1):
                y = [1, 1] + [t * ((k % 2) * 2 - 1) for k in range(n - 2)]
                def img(map_):
                    x = [map_[3][q] for q in range(d)]
                    for r in range(n):
                        x[map_[1][r]] += map_[2][r] * y[r]
                    return x
                self.assertEqual(img(a), img(b))

    def test_loop_words_need_two_double_crossings(self):
        # LEMMA: a nonempty word that ends in its start cell crosses at least two lattice axes twice, so a geodesic loop is at least
        # sqrt 2 (3h - u)... >= s sqrt 2 >= 1 long once a ball fits (s >= 1/sqrt 2). Exhaustive for d = 5 (all 109,601 words per cell) and,
        # with one crossing per axis, for d = 6..8 (none returns). The proof (BOUNDARY.md 4B) is the signed-rotation product argument.
        for c0 in (0, 1, 6):
            ws = words(5, c0, 2)
            loops = [(w, m) for (w, m, steps) in ws if steps and w[0] == c0]
            self.assertGreater(len(loops), 0)
            for w, m in loops:
                self.assertGreaterEqual(sum(1 for x in m if abs(x) == 2), 2)
        for d in (6, 7, 8):
            self.assertEqual([1 for (w, m, steps) in words(d, 2, 1) if steps and w[0] == 2], [])

    def test_reach_search_equals_full_enumeration(self):
        # pairs of random points, d = 5 with two crossings per axis (s from 0.75) and d = 6, 7 with one (s >= 1): the pruned search finds exactly
        # the realised words of the full enumeration whose length is within the reach, and its minimum is the same
        random.seed(5)
        for d, mm, lo, trials in ((5, 2, 75, 14), (6, 1, 100, 30), (7, 1, 100, 12)):
            full = {}
            n = d - 1
            for _ in range(trials):
                while True:
                    s = Fr(random.randint(lo, 160), 100)
                    h = s / 2
                    c1, c2 = random.randrange(2 * d), random.randrange(2 * d)
                    mk = lambda: [Fr(random.randint(-1000, 1000), 1000) * h for _ in range(n)]  # noqa: E731
                    cfg_ = N.Cfg(f"{d}cube", "ball", 2, s, [c1, c2], [mk(), mk()])
                    if sum((a - b) ** 2 for a, b in zip(cfg_.X[0], cfg_.X[1])) <= cfg_.D ** 2:
                        break
                H, D = cfg_.H, cfg_.D
                if c1 not in full:
                    full[c1] = words(d, c1, mm)
                for loop in (False, True):
                    X2 = cfg_.X[0] if loop else cfg_.X[1]
                    ends = set(N.cells_containing(d, X2, H))
                    brute = sorted(
                        L
                        for (w, m, steps) in full[c1]
                        if w[0] in ends and not (loop and not steps)
                        for L in [N.realised_sq(d, H, cfg_.U[0], X2, w, m, steps)]
                        if L is not None and L <= D * D
                    )
                    got = sorted(x for x in N.path_lengths_sq(d, H, c1, cfg_.U[0], X2, D * D, exclude_empty=loop) if x <= D * D)
                    self.assertEqual(got, brute)
                    # and with both bounds off
                    off = sorted(x for x in N.path_lengths_sq(d, H, c1, cfg_.U[0], X2, D * D, exclude_empty=loop, prune=False) if x <= D * D) if d == 5 else got
                    self.assertEqual(off, brute)

    def test_clearance_closed_form(self):
        # the nearest hinge of the development: a lattice point with two odd multiples of h among its coordinates, the other coordinates free;
        # brute force over the pairs of axes and the odd multiples -5h..5h
        random.seed(8)
        for d in (5, 6, 8, 10):
            n = d - 1
            for _ in range(60):
                h = Fr(random.randint(35, 120), 100)
                u = [Fr(random.randint(-1000, 1000), 1000) * h for _ in range(n)]
                best = min(
                    (u[a] - oa * h) ** 2 + (u[b] - ob * h) ** 2
                    for a in range(n)
                    for b in range(a + 1, n)
                    for oa in range(-5, 6, 2)
                    for ob in range(-5, 6, 2)
                )
                self.assertEqual(N.clearance_sq(u, h), best)

    def test_one_ball_bound_and_where_it_is_attained(self):
        # clearance <= h sqrt 2, so a ball fits only from s = 1/sqrt 2, and exactly there at the cell's centre
        for d in (5, 10):
            n = d - 1
            s = Fr(7071, 10000)  # just below 1/sqrt 2
            self.assertFalse(verdict(d, "ball", cfg(s, [(0, [0] * n)])).ok)
            s = Fr(70711, 100000)  # just above
            self.assertTrue(verdict(d, "ball", cfg(s, [(0, [0] * n)]), strict=True).ok)


class Packings(unittest.TestCase):
    def test_every_test_packing_certifies_strictly(self):
        data = load("ncube_packings.json")["packings"]
        self.assertEqual(len(data), 6 * 7)
        for p in data:
            cert = C.certify(p["id"], p["config"])
            self.assertTrue(cert["valid"], p["id"])
            self.assertTrue(verdict(int(p["polytope"][:-4]), p["piece"], p["config"], strict=False).ok)

    def test_shrunk_copy_fails(self):
        for p in load("ncube_packings.json")["packings"]:
            c = p["config"]
            f = 1 - 1e-6
            shrunk = {"s": c["s"] * f, "pieces": [{"cell": q["cell"], "u": [x * f for x in q["u"]]} for q in c["pieces"]]}
            self.assertFalse(verdict(int(p["polytope"][:-4]), p["piece"], shrunk).ok, p["id"])

    def test_cube_tilings(self):
        # P = 2d k^(d-1) cubes at s = k: k^(d-1) per cell on the grid of spacing 1: valid, touching, so not strictly valid; one more cube fails
        for d, k in ((5, 1), (5, 2), (6, 2), (7, 2), (8, 1), (10, 1)):
            n = d - 1
            pts = []
            for idx in range(k ** n):
                r = idx
                p = []
                for _ in range(n):
                    p.append(Fr(2 * (r % k) - (k - 1), 2))
                    r //= k
                pts.append(p)
            pieces = [(c, p) for c in range(2 * d) for p in pts]
            P = 2 * d * k ** n
            self.assertEqual(len(pieces), P)
            c = cfg(str(k), [(c_, [str(x) for x in p]) for c_, p in pieces])
            v = verdict(d, "cube", c)
            self.assertTrue(v.ok, (d, k, v.reasons[:2]))
            # touching exactly (with k = 1 each cell holds one cube: no pair of parts in a cell)
            self.assertEqual(v.margins["part_separation"], 0 if k > 1 else None)
            self.assertFalse(verdict(d, "cube", c, strict=True).ok)
            # the same at s = k (1 + 1e-9) with the centres scaled alike: strictly valid
            lam = Fr(10 ** 9 + 1, 10 ** 9)
            cs = cfg(str(k * lam), [(c_, [str(x * lam) for x in p]) for c_, p in pieces])
            self.assertTrue(verdict(d, "cube", cs, strict=True).ok)
            if P < 9999:
                dup = cfg(str(k), [(c_, [str(x) for x in p]) for c_, p in pieces] + [(0, ["0"] * n)])
                self.assertFalse(verdict(d, "cube", dup).ok)

    def test_ball_grids(self):
        # k^n balls per cell on the grid of spacing 1 at s = k, d = 5 (160 balls) and 6 (384): touching everywhere, valid, not strictly
        for d in (5, 6):
            n = d - 1
            k = 2
            pts = [[Fr(2 * ((idx >> a) & 1) - 1, 2) for a in range(n)] for idx in range(k ** n)]
            c = cfg("2", [(c_, [str(x) for x in p]) for c_ in range(2 * d) for p in pts])
            t0 = time.time()
            v = verdict(d, "ball", c)
            self.assertTrue(v.ok, v.reasons[:2])
            self.assertEqual(v.margins["distance_sq"], 1)
            self.assertFalse(verdict(d, "ball", c, strict=True).ok)
            print(f"  [{d}-cube, {len(c['pieces'])} balls: {time.time() - t0:.1f} s]", file=sys.stderr)

    def test_lower_bounds_against_the_engines(self):
        # the engine's bounds (written by the cross-check script) pass the exact test, and the next float up fails it
        data = load("ncube_bounds.json")["bounds"]
        self.assertGreater(len(data), 100)
        for r in data:
            checks = N.exact_bounds(f"{r['d']}cube", r["piece"], r["P"])
            x = Fr(r["value"])
            self.assertTrue(any(ok(x) for (_, _, ok) in checks), r)
            up = Fr(math.nextafter(r["value"], math.inf))
            self.assertFalse(any(ok(up) for (_, _, ok) in checks), r)


class Adversarial(unittest.TestCase):
    def test_ball_over_a_hinge_face(self):
        # two coordinates a from their walls: the hinge is at sqrt 2 a; a = 9/25 clears 1/2 (sqrt 2 * .36 = .509), a = 7/20 does not
        for d in range(5, 11):
            n = d - 1
            s = Fr(3)
            h = s / 2

            def ball(a, c=0):
                u = [Fr(0)] * n
                u[0], u[1] = h - a, h - a
                return cfg(s, [(c, u)])

            self.assertTrue(verdict(d, "ball", ball(Fr(9, 25)), strict=True).ok)
            bad = verdict(d, "ball", ball(Fr(7, 20)))
            self.assertFalse(bad.ok)
            self.assertIn("(d-3)-face", bad.reasons[0])
            self.assertFalse(verdict(d, "ball", cfg(s, [(0, [h, h] + [0] * (n - 2))])).ok)  # centred on the hinge itself
            # touching a wall alone (one slack 0) is fine while the other slacks keep 1/2 away
            self.assertTrue(verdict(d, "ball", cfg(s, [(0, [h] + [0] * (n - 1))]), strict=True).ok)

    def pair(self, d, eps):
        """Ball A in cell (0,+) at 1/2 from the ridge x1 = h, ball B in cell (1,+) at 1/2 - eps from the ridge x0 = h: unfolded 1 - eps apart (eps > 0 overlaps)."""
        n = d - 1
        s = Fr(3)
        h = s / 2
        a = [Fr(0)] * n  # cell 1 (x0 = h): local axes 1, 2, ...: u[0] is x1
        a[0] = h - Fr(1, 2)
        b = [Fr(0)] * n  # cell 3 (x1 = h): local axes 0, 2, ...: u[0] is x0
        b[0] = h - Fr(1, 2) + eps
        return cfg(s, [(1, a), (3, b)])

    def test_balls_across_a_ridge(self):
        eps = Fr(1, 10 ** 9)
        for d in range(5, 11):
            ok = verdict(d, "ball", self.pair(d, Fr(0)))
            self.assertTrue(ok.ok)  # touching: unfolded distance exactly 1
            self.assertEqual(ok.margins["distance_sq"], 1)
            self.assertFalse(verdict(d, "ball", self.pair(d, Fr(0)), strict=True).ok)
            self.assertFalse(verdict(d, "ball", self.pair(d, eps)).ok)  # 1e-9 too close
            self.assertTrue(verdict(d, "ball", self.pair(d, -eps), strict=True).ok)  # 1e-9 of room

    def test_a_point_on_a_ridge_is_seen_from_both_cells(self):
        # the same point listed in the two cells that share its ridge is an overlap; 0.999 along the ridge too; 1 along it is clear
        for d in range(5, 11):
            n = d - 1
            h = Fr(3, 2)

            def at(along):
                a = [Fr(0)] * n
                a[0] = h
                b = [Fr(0)] * n
                b[0], b[1] = h, along
                return cfg(3, [(1, a), (3, b)])

            self.assertFalse(verdict(d, "ball", at(Fr(0))).ok)
            self.assertFalse(verdict(d, "ball", at(Fr(999, 1000))).ok)
            self.assertTrue(verdict(d, "ball", at(Fr(1))).ok)
            self.assertFalse(verdict(d, "ball", at(Fr(1)), strict=True).ok)

    def test_balls_overlapping_around_a_hinge(self):
        # Round the hinge x0 = -h, x1 = h, x2 = h the three cells 0, (1,+), (2,+) fill 270 degrees. A in cell 0 at radius .6 from the hinge and
        # angle 85 degrees (nearly on the ridge to (1,+)), B in cell (2,+) at radius .6 and angle 185 degrees as seen through (1,+): the straight
        # path through (1,+) crosses two ridges and is 2 (.6) sin 50 degrees = .919 long, the other way, across the ridge of cells 0 and (2,+),
        # 170 degrees round and 2 (.6) sin 85 degrees = 1.195 long. Only the word that crosses two ridges sees the overlap.
        rad = math.radians
        for d in range(5, 11):
            n = d - 1
            h = Fr(3, 2)
            r = 0.6
            w = N.step_map(d, N.step_map(d, N.start_map(d, 0), 0, 1), 1, 1)
            u1 = [Fr(0)] * n
            u1[0], u1[1] = h - Fr(r * math.cos(rad(85))), h - Fr(r * math.sin(rad(85)))

            def target(phi):
                y = [Fr(0)] * n
                y[0], y[1] = h + Fr(r * math.cos(rad(phi - 180))), h + Fr(r * math.sin(rad(phi - 180)))
                x = [h * w[3][t] for t in range(d)]
                for k in range(n):
                    x[w[1][k]] += w[2][k] * y[k]
                return [x[a] for a in N.axes(d, w[0] >> 1)]

            bad = cfg(3, [(0, u1), (w[0], target(185))])
            v = verdict(d, "ball", bad)
            self.assertFalse(v.ok)
            self.assertEqual(len(v.reasons), 1)
            self.assertIn("balls 0 and 1 overlap", v.reasons[0])
            self.assertAlmostEqual(math.sqrt(v.margins["distance_sq"]), 2 * r * math.sin(rad(50)), places=3)
            # the unpruned search finds the same word, and the pair at 235 degrees (150 and 120 degrees round) is clear
            if d == 5:
                v2 = verdict(d, "ball", bad, prune=False)
                self.assertEqual(v2.margins["distance_sq"], v.margins["distance_sq"])
            fine = cfg(3, [(0, u1), (w[0], target(235))])
            self.assertTrue(verdict(d, "ball", fine, strict=True).ok)

    def test_cube_over_a_hinge(self):
        for d in range(5, 11):
            n = d - 1
            s = Fr(3)
            h = s / 2
            u = [Fr(0)] * n
            u[0], u[1] = h - Fr(1, 4), h - Fr(1, 4)  # crosses two faces: the hinge passes through it
            v = verdict(d, "cube", cfg(s, [(0, u)]))
            self.assertFalse(v.ok)
            self.assertIn("(d-3)-face", v.reasons[0])
            u[1] = Fr(0)  # one face only
            self.assertTrue(verdict(d, "cube", cfg(s, [(0, u)])).ok)
            # touching a face from inside counts as crossing it when strict: two touching axes are a hinge in the closed cube
            u2 = [Fr(0)] * n
            u2[0], u2[1] = h - Fr(1, 2), h - Fr(1, 2)
            self.assertTrue(verdict(d, "cube", cfg(s, [(0, u2)])).ok)
            self.assertFalse(verdict(d, "cube", cfg(s, [(0, u2)]), strict=True).ok)

    def test_cubes_across_a_ridge(self):
        # A in cell (0,+) crosses the face x1 = h by 1/4: its slab lies in cell (1,+) along x0 from h - 1/4 to h. B in cell (1,+) with its
        # x0-extent [h - 5/4 + delta, h - 1/4 + delta] meets the slab in a layer of thickness delta
        for d in range(5, 11):
            n = d - 1
            s = Fr(3)
            h = s / 2
            a = [Fr(0)] * n
            a[0] = h - Fr(1, 4)  # cell 1 (x0 = h): u[0] is x1; |u| + 1/2 - h = 1/4
            for delta, ok_ns, ok_s in ((Fr(0), True, False), (Fr(-1, 10 ** 9), True, True), (Fr(1, 10 ** 9), False, False)):
                b = [Fr(0)] * n
                b[0] = h - Fr(3, 4) + delta  # cell 3 (x1 = h): u[0] is x0
                c = cfg(s, [(1, a), (3, b)])
                self.assertEqual(verdict(d, "cube", c).ok, ok_ns, (d, delta))
                self.assertEqual(verdict(d, "cube", c, strict=True).ok, ok_s, (d, delta))

    def test_out_of_scope_and_malformed(self):
        d = 6
        good = cfg(3, [(0, [0] * 5)])
        bad = lambda p, c=good: C.check(p, c)  # noqa: E731
        with self.assertRaises(C.CertError):
            C.check(pid(d, "ball", 1), {"s": 3, "pieces": [{"cell": 0, "u": [0] * 5, "q": [1, 0, 0, 0]}]})
        with self.assertRaises(C.CertError):
            C.check(pid(d, "ball", 1), {"s": 3, "pieces": [{"cell": 12, "u": [0] * 5}]})
        with self.assertRaises(C.CertError):
            C.check(pid(d, "ball", 1), {"s": 3, "pieces": [{"cell": 0, "u": [0] * 3}]})
        with self.assertRaises(C.CertError):
            C.check(pid(d, "ball", 2), good)
        with self.assertRaises(C.CertError):
            C.check(pid(d, "ball", 1), {"s": float("inf"), "pieces": good["pieces"]})
        with self.assertRaises(C.CertError):
            C.check(pid(d, "ball", 1), {"s": 3, "pieces": [{"cell": 0, "u": [0, 0, 0, 0, "x"]}]})
        with self.assertRaises(C.CertError):
            C.check("boundary/11cube/ball/1", good)
        with self.assertRaises(C.CertError):
            C.check("boundary/6cube/sphere/1", good)
        # a centre outside its cell is invalid, not malformed
        out = verdict(d, "ball", cfg(3, [(0, [2, 0, 0, 0, 0])]))
        self.assertFalse(out.ok)
        self.assertIn("outside", out.reasons[0])
        # a cube outside the problem set's range of s: s < 1 crosses every axis
        self.assertFalse(verdict(d, "cube", cfg(Fr(9, 10), [(0, [0] * 5)])).ok)
        _ = bad

    def test_work_limit(self):
        # a search over its budget is refused with a plain error, never answered
        d = 8
        n = d - 1
        pieces = [(0, [Fr(1, 100) * i] * n) for i in range(2)]
        c = cfg(Fr(7071, 10000) * 2, pieces)
        with self.assertRaises(N.NCubeError):
            verdict(d, "ball", c, work=N.Work(50))


class CrossCheck(unittest.TestCase):
    def test_engine_verdicts(self):
        data = load("ncube_crosscheck.json")["cases"]
        self.assertGreater(len(data), 600)
        n_valid = 0
        borderline = 0
        t0 = time.time()
        for c in data:
            d = c["d"]
            v = verdict(d, c["piece"], {"s": c["s"], "pieces": c["pieces"]})
            m = v.margins
            near = False
            if c["piece"] == "ball":
                near = (m.get("clearance_sq") is not None and abs(m["clearance_sq"] - Fr(1, 4)) < Fr(1, 10 ** 11)) or (
                    m.get("distance_sq") is not None and abs(m["distance_sq"] - 1) < Fr(1, 10 ** 11)
                )
            else:
                near = (m.get("edge_slack") is not None and abs(m["edge_slack"]) < Fr(1, 10 ** 11)) or (m.get("part_separation") is not None and abs(m["part_separation"]) < Fr(1, 10 ** 11))
            if near:
                borderline += 1
            else:
                self.assertEqual(c["tol0"], v.ok, (c["source"], c["reason"], v.reasons[:1]))
            # the engine's tolerance 1e-9 only widens its verdict; an exactly invalid case far from valid stays invalid
            if v.ok:
                n_valid += 1
                self.assertTrue(c["tol9"], c["source"])
        self.assertGreater(n_valid, 250)
        print(f"  [{len(data)} engine verdicts, {n_valid} valid, {borderline} borderline: {time.time() - t0:.1f} s]", file=sys.stderr)


class Records(unittest.TestCase):
    def test_ids_and_front_end(self):
        self.assertTrue(C.ID_RE.match("boundary/10cube/cube/400"))
        self.assertTrue(C.ID_RE.match("boundary/5cube/ball/1"))
        self.assertIsNone(C.ID_RE.match("boundary/11cube/ball/1"))
        self.assertIsNone(C.ID_RE.match("boundary/5cube/ball/01"))
        p = load("ncube_packings.json")["packings"][1]
        out = C.certify(p["id"], p["config"])
        self.assertEqual(out["checker"], "python/boundary/ncube.py")
        # published bounds go through certify.lower_bounds
        kinds = [k for (k, _, _) in C.lower_bounds("ball", 20, "10cube")]
        self.assertEqual(kinds, ["volume", "one piece"])


if __name__ == "__main__":
    t0 = time.time()
    prog = unittest.main(exit=False, verbosity=1)
    print(f"test_ncube: {time.time() - t0:.1f} s", file=sys.stderr)
    sys.exit(0 if prog.result.wasSuccessful() else 1)
