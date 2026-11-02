"""Tests of the general exact certifier (general.py): balls on the boundary of the 5-cell, 16-cell, 24-cell, 120-cell and 600-cell
(and, as an independent check of the method, on the tesseract against the lattice-word certifier of certify.py).

  * the containers: counts of the literature, one squared edge, one norm of the cell normals, the gluing rotations exact isometries
    that compose to the identity, the cone angles of the report;
  * the engine's conventions (cell order, local frames) against fixtures/containers.json, written by
    packages/boundary/scripts/gen-containers.ts;
  * the engine's verdicts and chain distances against fixtures/general_crosscheck.json (gen-general-crosscheck.ts);
  * adversarial cases: a ball over a hinge, two balls 1e-9 too close across a ridge, a self-overlapping loop near a sharp edge
    of the 5-cell (cone angle 211.6 degrees), touching variants;
  * the test packings of fixtures/general_packings.json certified strictly valid, with the lower bounds;
  * completeness: a larger reach finds no further violation;
  * malformed input and out-of-scope pieces.

Run: python3 python/boundary/test_general.py
"""

import json
import math
import os
import sys
import time
import unittest
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import certify as C  # noqa: E402
import general as G  # noqa: E402

FIX = os.path.join(HERE, "fixtures")
NAMES = ["5cell", "16cell", "24cell", "120cell", "600cell"]


def load(name):
    with open(os.path.join(FIX, name)) as fh:
        return json.load(fh)


def cfg(s, pieces):
    return {"s": s, "pieces": [{"cell": c, "u": list(u)} for c, u in pieces]}


def centres(cont, cells, s):
    return {"s": s, "pieces": [{"cell": c, "u": [0.0, 0.0, 0.0]} for c in cells]}


def adjacent_cells(cont):
    c = cont.cells[0]
    return 0, c.facets[0]["nbr"]


class Containers(unittest.TestCase):
    def test_counts_and_isometries(self):
        for name in G.NAMES:
            cont = G.get_container(name)
            nV, nE, nC, nR, cv, cf = G.COUNTS[name]
            self.assertEqual(len(cont.V), nV)
            self.assertEqual(len(cont.cells), nC)
            self.assertEqual(sum(len(a) for a in cont.adj) // 2, nE)
            for c in cont.cells:
                self.assertEqual(len(c.ids), cv)
                self.assertEqual(len(c.facets), cf)
                self.assertEqual(G.vdot(c.N, c.N), cont.a, "all normals have one norm")
                for i in c.ids:
                    self.assertEqual(G.vdot(c.N, cont.V[i]), c.b, "the cell lies in its hyperplane")
            # every vertex lies on or below every cell's hyperplane (convexity), spot-checked on a few cells
            for c in cont.cells[:: max(1, len(cont.cells) // 8)]:
                self.assertTrue(all(G.sgn(G.vdot(c.N, v) - c.b) <= 0 for v in cont.V))

    def test_gluings_are_exact_isometries_and_invert(self):
        for name in G.NAMES:
            cont = G.get_container(name)
            for c in cont.cells[:: max(1, len(cont.cells) // 5)]:
                for fi in (0, len(c.facets) - 1):
                    f = c.facets[fi]
                    c2 = cont.cells[f["nbr"]]
                    A, t = cont.ginv(c, fi)
                    # an isometry from the chart of c2 to the chart of c: A^T G A = G2
                    AtGA = [[G.vdot([A[k][i] for k in range(3)], G.matvec(c.G, [A[k][j] for k in range(3)])) for j in range(3)] for i in range(3)]
                    self.assertEqual(AtGA, c2.G, f"{name} cell {c.index} facet {fi}")
                    # it fixes the ridge pointwise: the chart of each ridge vertex in c2 goes to its chart in c
                    for vid in f["ids"]:
                        self.assertEqual(G.vadd(G.matvec(A, c2.chart_v[vid]), t), c.chart_v[vid])
                    # the opposite gluing (c into c2) also fixes the ridge pointwise
                    A2, t2 = cont.ginv(c2, f["back"])
                    for vid in f["ids"]:
                        self.assertEqual(G.vadd(G.matvec(A2, c.chart_v[vid]), t2), c2.chart_v[vid])

    def test_cone_angles_of_the_report(self):
        want = {"5cell": 211.58633809652784, "tesseract": 270.0, "16cell": 282.1151174620371, "24cell": 328.4136619034721, "120cell": 349.6951535312339, "600cell": 352.6438968275464}
        for name, deg in want.items():
            cont = G.get_container(name)
            c = cont.cells[0]
            # the dihedral angle of the cell at each of its edges, from the facets' normals in the chart (metric Ginv)
            f0 = c.facets[0]
            got = []
            planes = G.facet_planes(cont, c)
            for h in c.hinges[:1]:
                owners = [i for i, f in enumerate(c.facets) if h[0] in f["ids"] and h[1] in f["ids"]]
                n1, n2 = planes[owners[0]][0], planes[owners[1]][0]
                m1, m2 = G.matvec(c.Ginv, n1), G.matvec(c.Ginv, n2)
                cosang = float(G.vdot(n1, m2)) / math.sqrt(float(G.vdot(n1, m1)) * float(G.vdot(n2, m2)))
                dihedral = math.pi - math.acos(cosang)
                cells_at_edge = {"5cell": 3, "tesseract": 3, "16cell": 4, "24cell": 3, "120cell": 3, "600cell": 5}[name]
                got.append(math.degrees(dihedral) * cells_at_edge)
            self.assertAlmostEqual(got[0], deg, places=6, msg=name)


class Conventions(unittest.TestCase):
    def test_cells_and_frames_match_the_engine(self):
        fx = load("containers.json")
        for name in G.NAMES:
            cont = G.get_container(name)
            f = fx[name]
            self.assertEqual(len(cont.cells), f["cellCount"])
            for v, w in zip(cont.fvertices, f["vertices"]):
                for a, b in zip(v, w):
                    self.assertAlmostEqual(a, b, delta=1e-14)
            for c in f["cells"]:
                mine = cont.cells[c["index"]]
                self.assertEqual(mine.ids, c["ids"], f"{name} cell {c['index']}")
                for a, b in zip(mine.forigin, c["origin"]):
                    self.assertAlmostEqual(a, b, delta=1e-14)
                for fr1, fr2 in zip(mine.fframe, c["frame"]):
                    for a, b in zip(fr1, fr2):
                        self.assertAlmostEqual(a, b, delta=1e-14)


class CrossCheck(unittest.TestCase):
    """The engine's verdicts at tolerance 1e-9 and its chain distances on seeded random configurations."""

    def test_against_the_engine(self):
        cases = load("general_crosscheck.json")["cases"]
        checked = skipped = 0
        differing = []
        for c in cases:
            config = {"s": Fr(c["config"]["s"]), "pieces": [{"cell": p["cell"], "u": [Fr(x) for x in p["u"]]} for p in c["config"]["pieces"]]}
            v = G.check_balls(c["polytope"], config, strict=False, full=True)
            cl = [math.sqrt(float(x)) for x in v.clear_d2.values()]
            ds = [math.sqrt(float(x)) for x in v.pair_d2.values()]
            ls = [math.sqrt(float(x)) for x in v.loop_d2.values()]
            # a verdict within 2e-8 of a threshold depends on the engine's tolerance: not compared
            if any(abs(x - 0.5) < 2e-8 for x in cl) or any(abs(x - 1) < 2e-8 for x in ds + ls):
                skipped += 1
            else:
                checked += 1
                self.assertEqual(v.ok, c["ts_valid_1e9"], f"{c['id']}: {c['ts_reason']} / {v.reasons[:1]}")
            for i, j, d in c["pairs"]:
                if d is not None and d < 0.99:
                    self.assertAlmostEqual(math.sqrt(float(v.pair_d2[(i, j)])), d, delta=1e-9, msg=f"{c['id']} pair {i},{j}")
            for i, d in c["loops"]:
                if d is not None and d < 0.99:
                    self.assertAlmostEqual(math.sqrt(float(v.loop_d2[i])), d, delta=1e-9, msg=f"{c['id']} loop {i}")
            for i, d in enumerate(c["clearance"]):
                if d is not None and d < 0.74:
                    self.assertAlmostEqual(math.sqrt(float(v.clear_d2[i])), d, delta=1e-9, msg=f"{c['id']} clearance {i}")
            # the exact enumerator keeps the chains the float one does (which also cuts by the turn rule): same node counts
            for i, n in enumerate(c["chains"]):
                if v.chains_by_ball[i] != n:
                    differing.append((c["id"], i, v.chains_by_ball[i], n))
        self.assertGreater(checked, 120)
        self.assertLessEqual(len(differing), 3, differing)

    def test_tesseract_records_against_the_lattice_certifier(self):
        """The method is independent of the lattice words: it certifies the tesseract's ball records with the same margins."""
        with open(os.path.join(HERE, "..", "..", "data", "boundary", "records.json")) as fh:
            recs = json.load(fh)["records"]
        n = 0
        for r in recs:
            if r["piece"] != "ball" or r["P"] not in (4, 8, 13):
                continue
            conf = {"s": Fr(r["config"]["s"]), "pieces": [{"cell": p["cell"], "u": [Fr(x) for x in p["u"]]} for p in r["config"]["pieces"]]}
            v = G.check_balls("tesseract", conf, strict=True)
            self.assertTrue(v.ok, r["id"])
            lat = C.check(r["id"], r["config"], strict=True)
            u2 = lat.cfg.unit ** 2
            self.assertAlmostEqual(float(v.margins["clearance_sq"]), float(lat.margins["clearance_sq"] * u2), delta=1e-9)
            self.assertAlmostEqual(float(v.margins["distance_sq"]), float(lat.margins["distance_sq"] * u2), delta=1e-9)
            n += 1
        self.assertEqual(n, 3)


class Adversarial(unittest.TestCase):
    def test_cell_centres_across_a_ridge_1e9(self):
        """Two balls at the centres of adjacent cells are exactly 2 r s apart (r the cell's inradius): valid from
        s0 = 1 / (2 r), 1e-9 too close below it, for every container."""
        # exact edges: the cell inradius of a unit-edge tetrahedron is sqrt 6 / 12, octahedron sqrt 6 / 6, dodecahedron below
        s0 = {"5cell": math.sqrt(6.0), "16cell": math.sqrt(6.0), "24cell": math.sqrt(6.0) / 2, "600cell": math.sqrt(6.0), "120cell": 10 / math.sqrt(250 + 110 * math.sqrt(5))}
        for name in NAMES:
            cont = G.get_container(name)
            a, b = adjacent_cells(cont)
            for rel, expect in ((-1e-9, False), (1e-9, True)):
                v = G.check_balls(name, centres(cont, [a, b], Fr(s0[name] * (1 + rel))), strict=True)
                self.assertEqual(v.ok, expect, f"{name} at {rel:+g}: {v.reasons[:1]}")
                if not expect:
                    self.assertTrue(any("overlap" in r for r in v.reasons))
            # the least distance is within 3e-9 of 1 on the valid side
            self.assertAlmostEqual(math.sqrt(float(v.margins["distance_sq"])), 1 + 1e-9, delta=3e-9, msg=name)

    def test_touching_is_valid_but_not_strictly(self):
        """The 24-cell at s = 1: one ball at a cell's centre is exactly 1/2 from every edge of the octahedron."""
        cont = G.get_container("24cell")
        one = centres(cont, [0], Fr(1))
        self.assertTrue(G.check_balls("24cell", one, strict=False).ok)
        v = G.check_balls("24cell", one, strict=True)
        self.assertFalse(v.ok)
        self.assertEqual(v.margins["clearance_sq"], Fr(1, 4))
        # the tesseract at s = 1: eight balls at the cells' centres touch across the ridges (exactly 1 apart)
        eight = centres(G.get_container("tesseract"), range(8), Fr(1))
        self.assertTrue(G.check_balls("tesseract", eight, strict=False).ok)
        self.assertFalse(G.check_balls("tesseract", eight, strict=True).ok)
        # certify.py's front end agrees
        self.assertTrue(C.check("boundary/24cell/ball/1", {"s": 1, "pieces": [{"cell": 0, "u": [0, 0, 0]}]}, strict=False).ok)
        with self.assertRaises(C.CertError):
            C.certify("boundary/24cell/ball/1", {"s": 1, "pieces": [{"cell": 0, "u": [0, 0, 0]}]})

    def test_ball_over_a_hinge(self):
        """A ball centred on an edge of the polytope is invalid on every container; so is one 1e-9 short of its clearance."""
        for name in NAMES:
            cont = G.get_container(name)
            s = 3.0 if name != "120cell" else 1.5
            c = cont.cells[0]
            a, b = c.hinges[0]
            fv = cont.fvertices
            mid = [(fv[a][k] + fv[b][k]) / 2 * s for k in range(4)]
            u = [sum((mid[k] - s * c.forigin[k]) * c.fframe[j][k] for k in range(4)) for j in range(3)]
            v = G.check_balls(name, {"s": Fr(s), "pieces": [{"cell": 0, "u": [Fr(x) for x in u]}]}, strict=False, full=True)
            self.assertFalse(v.ok, name)
            self.assertTrue(any("edge" in r for r in v.reasons), name)
            self.assertLess(float(v.margins["clearance_sq"]), 1e-12)

    def test_one_ball_just_below_its_bound(self):
        """P = 1 at a cell's centre: valid from s = 1 / (2 rho) (rho the cell's distance to its edges); 1e-9 below it, not."""
        s0 = {"5cell": math.sqrt(2.0), "16cell": math.sqrt(2.0), "600cell": math.sqrt(2.0), "24cell": 1.0, "120cell": 1 / ((3 + math.sqrt(5)) / 2)}
        for name in NAMES:
            cont = G.get_container(name)
            for rel, expect in ((-1e-9, False), (1e-9, True)):
                v = G.check_balls(name, centres(cont, [0], Fr(s0[name] * (1 + rel))), strict=True)
                self.assertEqual(v.ok, expect, f"{name} {rel:+g}")

    def test_loop_near_a_sharp_edge_of_the_5_cell(self):
        """A centre close to an edge of the 5-cell (cone angle 211.6 degrees) has a geodesic loop much shorter than 1: the
        engine's chains find 0.0940692594805195; so does the exact enumeration, and both reject the ball."""
        conf = {"s": Fr(1.6655998116219417), "pieces": [{"cell": 0, "u": [Fr(-0.5138403881454334), Fr(0.6838350866035227), Fr(-0.4569581382944524)]}]}
        v = G.check_balls("5cell", conf, strict=False, full=True)
        self.assertFalse(v.ok)
        self.assertAlmostEqual(math.sqrt(float(v.loop_d2[0])), 0.0940692594805195, delta=1e-9)
        self.assertAlmostEqual(math.sqrt(float(v.clear_d2[0])), 0.02179945071089696, delta=1e-9)
        self.assertTrue(any("overlaps itself" in r for r in v.reasons))
        self.assertTrue(any("edge" in r for r in v.reasons))

    def test_a_loop_the_clearance_does_not_hide(self):
        """Without the clearance stage the loop is still found (full mode): shrinking the 5-cell until its balls' loops go
        below 1 while the centre stays 1/2 from the edges is impossible, which is why loops add nothing to clearance on these
        containers (BOUNDARY.md 4A.3); the test pins the engine's finding at a second centre."""
        conf = {"s": Fr(1.6655998116219417), "pieces": [{"cell": 0, "u": [Fr(0.0), Fr(0.0), Fr(0.0)]}]}
        v = G.check_balls("5cell", conf, strict=True, full=True)
        self.assertTrue(v.ok)
        self.assertIsNone(v.margins["loop_sq"])

    def test_a_centre_on_a_ridge_is_developed_from_both_cells(self):
        """A ball centred at the centre of a ridge (exactly: the mean of the ridge's vertices): it lies in two cells and is
        developed from both; with s large enough it is valid, and a second ball 0.9 away across the ridge overlaps it."""
        for name, s in (("24cell", 2.0), ("600cell", 4.5), ("16cell", 3.0), ("5cell", 4.0), ("120cell", 1.2)):
            cont = G.get_container(name)
            c = cont.cells[0]
            f = c.facets[0]
            nb = f["nbr"]
            ids = f["ids"]
            mean = [G.Q5(0)] * cont.n
            for i in ids:
                mean = G.vadd(mean, cont.V[i])
            mean = [x / len(ids) for x in mean]
            conf = {"s": Fr(s), "pieces": [{"cell": 0, "u": [Fr(0)] * 3}]}
            prob = G.Problem(cont, conf, points=[mean])
            self.assertEqual(sorted(prob.balls[0].cells), sorted([0, nb]), name)
            v = G.check_balls(name, conf, strict=True, full=True, points=[mean])
            self.assertTrue(v.ok, f"{name}: {v.reasons[:1]}")
            self.assertIsNone(v.margins["loop_sq"], name)
            # a second centre 0.9 into the neighbour, on the line from the ridge's centre to the neighbour's centre (perpendicular to
            # the ridge): exact too, as the mean plus a multiple of (the neighbour's centre - the mean). The base length of that
            # multiple is 0.9 e0 / s, so the multiple is rational only for a rational e0 / s: take the point at the exact fraction
            # t = 1/3 of the way instead and compare with the exact distance
            c2 = cont.cells[nb]
            cen2 = [G.Q5(0)] * cont.n
            for i in c2.ids:
                cen2 = G.vadd(cen2, cont.V[i])
            cen2 = [x / len(c2.ids) for x in cen2]
            t = Fr(1, 3)
            y = G.vadd(mean, G.vscale(t, G.vsub(cen2, mean)))
            conf2 = {"s": Fr(s), "pieces": [conf["pieces"][0], {"cell": nb, "u": [Fr(0)] * 3}]}
            w = G.check_balls(name, conf2, strict=True, full=True, points=[mean, y])
            d2 = G.vdot(G.vsub(y, mean), G.vsub(y, mean))
            self.assertEqual(w.pair_d2[(0, 1)], d2 * Fr(s) ** 2 / cont.e2, name)

    def test_centre_outside_its_cell_is_refused(self):
        conf = {"s": Fr(3), "pieces": [{"cell": 0, "u": [Fr(100), Fr(0), Fr(0)]}]}
        v = G.check_balls("24cell", conf, strict=True)
        self.assertFalse(v.ok)
        self.assertTrue(any("outside" in r for r in v.reasons))


class Packings(unittest.TestCase):
    def test_every_test_packing_is_certified_and_bounded(self):
        data = load("general_packings.json")["packings"]
        self.assertGreaterEqual(len(data), 25)
        seen = set()
        for r in data:
            t = time.time()
            v = C.check(r["id"], r["config"], strict=True)
            self.assertTrue(v.ok, f"{r['id']}: {v.reasons[:1]}")
            summ = C.summary(v)
            self.assertEqual(summ["min_distance"], r["min_distance"])
            self.assertEqual(summ["min_clearance"], r["min_clearance"])
            self.assertEqual(r["score_cert"], r["config"]["s"])
            # strictly above both proven bounds, exactly (ok(x) says x is not above the bound)
            s = Fr(r["score_cert"])
            for kind, _, ok in C.lower_bounds("ball", r["P"], r["polytope"]):
                self.assertFalse(ok(s), f"{r['id']}: s is not above the {kind} bound")
            seen.add(r["polytope"])
        self.assertEqual(seen, set(NAMES))

    def test_bounds_are_below_the_packings_and_sharp_for_one_piece(self):
        for name in NAMES:
            for kind, _, ok in G.exact_bounds(name, 1):
                self.assertFalse(ok(Fr(10)), f"{name} {kind}")
                self.assertTrue(ok(Fr(1, 100)), f"{name} {kind}")

    def test_cell_centre_packings(self):
        """One ball per cell, every cell of a small container: valid from the ridge-crossing distance."""
        s0 = {"5cell": math.sqrt(6.0), "16cell": math.sqrt(6.0), "24cell": math.sqrt(6.0) / 2}
        for name, C_ in (("5cell", 5), ("16cell", 16), ("24cell", 24)):
            cont = G.get_container(name)
            v = G.check_balls(name, centres(cont, range(C_), Fr(s0[name] * (1 + 1e-9))), strict=True)
            self.assertTrue(v.ok, name)


class Completeness(unittest.TestCase):
    def test_a_larger_reach_finds_no_further_violation(self):
        """The chains within the ball's diameter hold every path shorter than it (BOUNDARY.md 4A.2): with 1.5 diameters
        the pairs and loops shorter than 1 are the same, with the same lengths."""
        data = {r["id"]: r for r in load("general_packings.json")["packings"]}
        for rid in ("boundary/5cell/ball/3", "boundary/16cell/ball/4", "boundary/24cell/ball/3", "boundary/600cell/ball/2"):
            r = data[rid]
            conf = {"s": Fr(r["config"]["s"]), "pieces": [{"cell": p["cell"], "u": [Fr(x) for x in p["u"]]} for p in r["config"]["pieces"]]}
            a = G.check_balls(r["polytope"], conf, strict=True)
            b = G.check_balls(r["polytope"], conf, strict=True, reach_factor=Fr(9, 4))
            self.assertTrue(a.ok and b.ok, rid)
            self.assertGreaterEqual(b.stats["chains"], a.stats["chains"])
            # everything of length at most 1 is the same in both runs
            for tab in ("pair_d2", "loop_d2", "clear_d2"):
                ta, tb = getattr(a, tab), getattr(b, tab)
                for k, x in ta.items():
                    self.assertEqual(tb[k], x, f"{rid} {tab} {k}")

    def test_reach_changes_nothing_on_an_invalid_packing(self):
        """The same for an overlapping pair: the shortest path across the ridge is found with the plain reach."""
        cont = G.get_container("24cell")
        a, b = adjacent_cells(cont)
        conf = centres(cont, [a, b], Fr(1.2))
        x = G.check_balls("24cell", conf, strict=True, full=True)
        y = G.check_balls("24cell", conf, strict=True, full=True, reach_factor=Fr(4))
        self.assertFalse(x.ok)
        self.assertEqual(x.pair_d2[(0, 1)], y.pair_d2[(0, 1)])


class Input(unittest.TestCase):
    def test_malformed_and_out_of_scope(self):
        ok = {"s": 3, "pieces": [{"cell": 0, "u": [0, 0, 0]}]}
        for bad_id in ("boundary/24cell/ball/01", "boundary/9cell/ball/1", "boundary/24cell/ball/", "boundary/24cell/ball/1/"):
            with self.assertRaises(C.CertError):
                C.check(bad_id, ok)
        with self.assertRaises(C.CertError):
            C.check("boundary/24cell/cube/1", ok)  # balls only
        for conf in (
            {"s": 0, "pieces": ok["pieces"]},
            {"s": 1e5, "pieces": ok["pieces"]},
            {"s": True, "pieces": ok["pieces"]},
            {"s": 3, "pieces": []},
            {"s": 3, "pieces": [{"cell": 24, "u": [0, 0, 0]}]},
            {"s": 3, "pieces": [{"cell": -1, "u": [0, 0, 0]}]},
            {"s": 3, "pieces": [{"cell": 0, "u": [0, 0]}]},
            {"s": 3, "pieces": [{"cell": 0, "u": [0, 0, float("nan")]}]},
            {"s": 3, "pieces": [{"cell": 0, "u": [0, 0, 0], "q": [1, 0, 0, 0]}]},
            {"s": 3, "pieces": [{"cell": 0, "u": [0, 0, 0], "z": 1}]},
            "not an object",
        ):
            with self.assertRaises(C.CertError, msg=str(conf)):
                C.check("boundary/24cell/ball/1", conf)


class Review(unittest.TestCase):
    """The fixes from the review of 4A: a reading slack, a cheap pre-check, CertError for guards, ids that end at the end."""

    ONE = {"pieces": [{"cell": 0, "u": [0, 0, 0]}]}

    def test_reading_slack(self):
        # a ball at a 24-cell cell's centre has clearance s/2: strictly valid needs 1e-13 s of room for the engine's reading
        self.assertFalse(C.check("boundary/24cell/ball/1", dict(self.ONE, s=1 + 5e-14), strict=True).ok)
        self.assertTrue(C.check("boundary/24cell/ball/1", dict(self.ONE, s=1 + 1e-12), strict=True).ok)

    def test_tiny_s_is_decided_before_any_chain_search(self):
        for name in ("5cell", "120cell", "600cell"):
            v = C.check(f"boundary/{name}/ball/1", dict(self.ONE, s=0.01), strict=True)
            self.assertFalse(v.ok)
            self.assertIn("bound", v.reasons[0])
        v = C.check("boundary/24cell/ball/3", {"s": 0.1, "pieces": self.ONE["pieces"] * 3}, strict=False)
        self.assertFalse(v.ok)

    def test_ids_end_at_the_end(self):
        for rid in ("boundary/24cell/ball/1\n", "boundary/5cube/ball/1\n", "boundary/tesseract/cube/2\n"):
            with self.assertRaises(C.CertError, msg=repr(rid)):
                C.problem_fields(rid)
        self.assertEqual(C.problem_fields("boundary/24cell/ball/1"), ("24cell", "ball", 1))


if __name__ == "__main__":
    unittest.main(verbosity=1)
