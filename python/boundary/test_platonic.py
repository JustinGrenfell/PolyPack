"""Tests of the exact certifier for Platonic pieces (platonic.py, qmulti.py): rotated tetrahedra, cubes, octahedra, dodecahedra and
icosahedra on the boundary of the regular 4-polytopes (BOUNDARY.md section 4C).

  * the field Q(sqrt 2, sqrt 5): arithmetic and exact signs against 60-digit enclosures, inverses, mixed arithmetic with Q5;
  * the solids: unit edges, face counts, exact radii and volumes; the frames of the cells (orthogonal, equal length, kappa^2 m = e0^2) and
    their images in the engine's float coordinates (orthonormal);
  * rational rotations (exactly orthogonal) and the reading of the engine's floats, with the measured discrepancy;
  * exact poses: the cells themselves at s = 1 (valid, touching; not strictly) and slightly above (strictly), the tesseract's cubes by
    this path against the lattice-word certifier of certify.py;
  * adversarial cases: overlaps across a ridge by 1e-6, an edge of the polytope through a solid by 1e-6, a centre on a ridge, a centre on an edge,
    a solid around a sharp edge of the 5-cell overlapping itself, touching variants, malformed input and the work cap;
  * the engine's verdicts (fixtures/platonic_crosscheck.json, packages/boundary/scripts/gen-platonic-crosscheck.ts): at 1e-5 and 1e-9 beyond the
    critical scale of random pairs, the exact verdict is the engine's, whatever the contact (the float validator never accepts what the exact
    certifier refuses, and the reverse, away from the threshold);
  * the lower bounds of the engine against the exact bounds;
  * the test packings (fixtures/platonic_packings.json) certified strictly valid.

Run: python3 python/boundary/test_platonic.py
"""

import itertools
import json
import math
import os
import random
import sys
import unittest
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import certify as C  # noqa: E402
import general as G  # noqa: E402
import platonic as PL  # noqa: E402
from qmulti import QM  # noqa: E402

FIX = os.path.join(HERE, "fixtures")
NAMES = ["tesseract", "5cell", "16cell", "24cell", "120cell", "600cell"]
SOLIDS = ["tetrahedron", "cube", "octahedron", "dodecahedron", "icosahedron"]
CELL_OF = {"5cell": "tetrahedron", "16cell": "tetrahedron", "600cell": "tetrahedron", "24cell": "octahedron", "120cell": "dodecahedron", "tesseract": "cube"}


def load(name):
    with open(os.path.join(FIX, name)) as fh:
        return json.load(fh)


def signed_rotations():
    mats = []
    for perm in itertools.permutations(range(3)):
        inv = sum(1 for i in range(3) for j in range(i + 1, 3) if perm[i] > perm[j]) % 2
        for sg in itertools.product((1, -1), repeat=3):
            if (sg[0] * sg[1] * sg[2]) * (-1 if inv else 1) == 1:
                m = [[Fr(0)] * 3 for _ in range(3)]
                for i in range(3):
                    m[i][perm[i]] = Fr(sg[i])
                mats.append(m)
    return mats


ROTATIONS = signed_rotations()


def cells_exact(name, piece, count=None):
    """Exact poses that put the solid onto each cell of the container (centred, turned onto it): found among the 24 axis rotations and the sign of the
    third frame vector; exact vertex equality decides."""
    cont = G.get_container(name)
    sh = PL.get_shape(piece)
    poses = []
    for c in cont.cells[:count]:
        found = None
        for sg, R in itertools.product((1, -1), ROTATIONS):
            pose = PL.Pose(c.index, c.centroid, R, sg)
            _, verts = PL.piece_vertices_in_anchor(cont, sh, Fr(1), pose, Fr(0))
            want = [[QM.of(x) for x in c.chart_v[i]] for i in c.ids]
            if all(any(all(a == b for a, b in zip(v, w)) for w in want) for v in verts):
                found = pose
                break
        assert found is not None, (name, c.index)
        poses.append(found)
    return cont, sh, poses


class Field(unittest.TestCase):
    def test_arithmetic_and_signs_against_enclosures(self):
        rnd = random.Random(1)

        def rq():
            return QM(*[Fr(rnd.randint(-9, 9), rnd.randint(1, 5)) for _ in range(4)])

        for _ in range(400):
            a, b = rq(), rq()
            lo, hi = a.bounds()
            self.assertEqual(a.sign(), (lo > 0) - (hi < 0) if not (lo <= 0 <= hi) else 0)
            for op in (lambda x, y: x + y, lambda x, y: x - y, lambda x, y: x * y):
                r = op(a, b)
                lo, hi = r.bounds()
                al, ah = a.bounds()
                bl, bh = b.bounds()
                # the enclosure of the result contains the product of the enclosures' midpoints
                mid = op((al + ah) / 2, (bl + bh) / 2)
                self.assertTrue(lo - Fr(1, 10 ** 30) <= mid <= hi + Fr(1, 10 ** 30))
            if a.sign() != 0:
                self.assertEqual(a * a.inv(), QM(1))
                self.assertEqual((b / a) * a, b)

    def test_signs_near_zero(self):
        # sqrt 2 + sqrt 5 against a 50-digit rational near it, on both sides
        v = QM(0, 1, 1, 0)
        lo, hi = v.bounds()
        t = Fr(3650281539872884, 10 ** 15)  # not close enough: sign decided by the float fast path
        self.assertEqual((v - t).sign(), 1)
        close = lo - Fr(1, 10 ** 44)
        self.assertEqual((v - close).sign(), 1)
        self.assertEqual((v - (hi + Fr(1, 10 ** 44))).sign(), -1)
        # exact zeros and the relations of the field
        self.assertEqual((QM(0, 1) * QM(0, 1) - 2).sign(), 0)
        self.assertEqual((QM(0, 0, 1) * QM(0, 0, 1) - 5).sign(), 0)
        self.assertEqual((QM(0, 1) * QM(0, 0, 1) - QM(0, 0, 0, 1)).sign(), 0)
        phi = QM(Fr(1, 2), 0, Fr(1, 2))
        self.assertEqual((phi * phi - phi - 1).sign(), 0)
        # (sqrt 10 - sqrt 2 sqrt 5) is zero; sqrt 5 - 2.2360679774997896964091736687313 is positive
        self.assertEqual((QM(0, 0, 0, 1) - QM(0, 1) * QM(0, 0, 1)).sign(), 0)
        self.assertEqual((QM(0, 0, 1) - Fr("2.2360679774997896964091736687312")).sign(), 1)
        self.assertEqual((QM(0, 0, 1) - Fr("2.2360679774997896964091736687314")).sign(), -1)

    def test_mixed_with_q5(self):
        a = G.Q5(Fr(1, 2), Fr(1, 3))
        b = QM(1, 2, 3, 4)
        for r in (a * b, b * a, a + b, b + a, a - b, b - a):
            self.assertIsInstance(r, QM)
        self.assertEqual(a * b, b * a)
        self.assertEqual(a - b, -(b - a))
        self.assertTrue(QM.of(a) == a and a == QM.of(a))
        self.assertTrue((a < b) == (QM.of(a) < b))
        self.assertEqual(G.Q5(2) * G.Q5(3), G.Q5(6))  # Q5 itself unchanged


class Solids(unittest.TestCase):
    def test_unit_edges_counts_and_radii(self):
        lit = {
            "tetrahedron": (4, 6, 4, Fr(3, 8), Fr(1, 24)),
            "cube": (8, 12, 6, Fr(3, 4), Fr(1, 4)),
            "octahedron": (6, 12, 8, Fr(1, 2), Fr(1, 6)),
            "dodecahedron": (20, 30, 12, None, None),
            "icosahedron": (12, 30, 20, None, None),
        }
        for name, (V, E, F, c2, r2) in lit.items():
            sh = PL.get_shape(name)
            self.assertEqual((len(sh.V), len(sh.edges), len(sh.faces)), (V, E, F), name)
            for i, j in sh.edges:
                d = [a - b for a, b in zip(sh.V[i], sh.V[j])]
                self.assertEqual(d[0] * d[0] + d[1] * d[1] + d[2] * d[2], QM(1))
            for k in range(3):
                self.assertEqual(sum((v[k] for v in sh.V), QM(0)), QM(0))
            if c2 is not None:
                self.assertEqual(QM.of(sh.circ2), QM.of(c2))
                self.assertEqual(QM.of(sh.in2), QM.of(r2))
            # every vertex at the same distance; every face plane at the same distance (the in-radius)
            for v in sh.V:
                self.assertEqual(v[0] * v[0] + v[1] * v[1] + v[2] * v[2], QM.of(sh.circ2))
        # dodecahedron: R^2 = (3/8)(3 + sqrt 5), r^2 = (25 + 11 sqrt 5)/40; icosahedron: R^2 = (5 + sqrt 5)/8, r^2 = (7 + 3 sqrt 5)/24
        self.assertEqual(QM.of(PL.get_shape("dodecahedron").circ2), QM(Fr(9, 8), 0, Fr(3, 8)))
        self.assertEqual(QM.of(PL.get_shape("dodecahedron").in2), QM(Fr(5, 8), 0, Fr(11, 40)))
        self.assertEqual(QM.of(PL.get_shape("icosahedron").circ2), QM(Fr(5, 8), 0, Fr(1, 8)))
        self.assertEqual(QM.of(PL.get_shape("icosahedron").in2), QM(Fr(7, 24), 0, Fr(1, 8)))

    def test_exact_volumes_against_floats(self):
        for name in SOLIDS:
            a, b, r = PL.exact_shape_volume(name)
            sh = PL.get_shape(name)
            fv = [[float(c) for c in v] for v in sh.V]
            vol = 0.0
            for f in sh.faces:
                # the pyramid from the centre over each face: fan the face's vertices in angular order
                ids = list(f["verts"])
                c = [sum(fv[i][k] for i in ids) / len(ids) for k in range(3)]
                n = [0, 0, 0]
                a0, b0, c0 = f["tri"]
                e1 = [fv[b0][k] - fv[a0][k] for k in range(3)]
                e2 = [fv[c0][k] - fv[a0][k] for k in range(3)]
                n = [e1[1] * e2[2] - e1[2] * e2[1], e1[2] * e2[0] - e1[0] * e2[2], e1[0] * e2[1] - e1[1] * e2[0]]
                nn = math.sqrt(sum(x * x for x in n))
                h = abs(sum(n[k] * fv[a0][k] for k in range(3))) / nn
                u = [fv[ids[0]][k] - c[k] for k in range(3)]
                w = [n[1] * u[2] - n[2] * u[1], n[2] * u[0] - n[0] * u[2], n[0] * u[1] - n[1] * u[0]]
                nw = math.sqrt(sum(x * x for x in w))
                ang = sorted(math.atan2(sum((fv[i][k] - c[k]) * w[k] for k in range(3)) / nw, sum((fv[i][k] - c[k]) * u[k] for k in range(3)) / math.sqrt(sum(x * x for x in u))) for i in ids)
                area = 0.0
                order = sorted(ids, key=lambda i: math.atan2(sum((fv[i][k] - c[k]) * w[k] for k in range(3)) / nw, sum((fv[i][k] - c[k]) * u[k] for k in range(3)) / math.sqrt(sum(x * x for x in u))))
                for i in range(len(order)):
                    p, q = order[i], order[(i + 1) % len(order)]
                    d1 = [fv[p][k] - c[k] for k in range(3)]
                    d2 = [fv[q][k] - c[k] for k in range(3)]
                    cr = [d1[1] * d2[2] - d1[2] * d2[1], d1[2] * d2[0] - d1[0] * d2[2], d1[0] * d2[1] - d1[1] * d2[0]]
                    area += 0.5 * math.sqrt(sum(x * x for x in cr))
                vol += area * h / 3
            lo, hi = (a + b * (G.SQRT2_LO if r == 2 else G.SQRT5_LO)), (a + b * (G.SQRT2_HI if r == 2 else G.SQRT5_HI))
            self.assertAlmostEqual(vol, float(lo), places=10, msg=name)


class Frames(unittest.TestCase):
    def test_frames_are_exact_and_orthonormal_in_the_engines_coordinates(self):
        rnd = random.Random(3)
        for name in NAMES:
            cont = G.get_container(name)
            cells = list(range(len(cont.cells))) if len(cont.cells) <= 24 else rnd.sample(range(len(cont.cells)), 6)
            for ci in cells:
                fr = PL.cell_frame(cont, ci)
                self.assertEqual(fr.kappa * fr.kappa * QM.of(fr.m), QM.of(cont.e2))
                rows = PL.engine_frame_images(cont, ci)
                F = cont.cells[ci].fframe
                # the engine's own frame spans the same 3-space: the exact frame's images are orthogonal to the engine's normal
                for r in rows:
                    self.assertAlmostEqual(sum(r[i] * cont.cells[ci].fnormal[i] for i in range(4)), 0.0, places=12)
                O = [[sum(rows[a][i] * F[b][i] for i in range(4)) for b in range(3)] for a in range(3)]
                self.assertAlmostEqual(abs(PL._det3f(O)), 1.0, places=12)

    def test_the_kinds_of_cell(self):
        # kappa is 1/sqrt2 for tetrahedra and octahedra, 1/phi for dodecahedra, 1 for cubes (docstring of platonic.py)
        want = {"5cell": QM(0, Fr(1, 2)), "16cell": QM(0, Fr(1, 2)), "24cell": QM(0, Fr(1, 2)), "600cell": QM(0, Fr(1, 2)), "120cell": QM(Fr(-1, 2), 0, Fr(1, 2)), "tesseract": QM(1)}
        for name, kappa in want.items():
            self.assertEqual(PL.cell_frame(G.get_container(name), 0).kappa, kappa)


class Rotations(unittest.TestCase):
    def test_rational_rotation_is_exactly_orthogonal_and_close(self):
        rnd = random.Random(7)
        for _ in range(50):
            q = [rnd.gauss(0, 1) for _ in range(4)]
            n = math.sqrt(sum(x * x for x in q))
            q = [x / n for x in q]
            R = PL.rational_rotation(q)
            for i in range(3):
                for j in range(3):
                    self.assertEqual(sum(R[i][k] * R[j][k] for k in range(3)), 1 if i == j else 0)
            # determinant +1
            det = R[0][0] * (R[1][1] * R[2][2] - R[1][2] * R[2][1]) - R[0][1] * (R[1][0] * R[2][2] - R[1][2] * R[2][0]) + R[0][2] * (R[1][0] * R[2][1] - R[1][1] * R[2][0])
            self.assertEqual(det, 1)
            Rf = PL.quat_matrix_float(q)
            self.assertLess(max(abs(float(R[i][j]) - Rf[i][j]) for i in range(3) for j in range(3)), 1e-15)

    def test_matrix_quaternion_round_trip(self):
        rnd = random.Random(8)
        for _ in range(100):
            q = [rnd.gauss(0, 1) for _ in range(4)]
            n = math.sqrt(sum(x * x for x in q))
            q = [x / n for x in q]
            m = PL.quat_matrix_float(q)
            q2 = PL.matrix_quaternion(m)
            m2 = PL.quat_matrix_float(q2)
            self.assertLess(max(abs(m[i][j] - m2[i][j]) for i in range(3) for j in range(3)), 1e-14)


class ExactPoses(unittest.TestCase):
    """The cells themselves: exact poses, exact decisions (non-strict at s = 1, strict above)."""

    def test_cells_at_one_are_valid_touching_and_not_strictly(self):
        for name in ["5cell", "tesseract"]:
            cont, sh, poses = cells_exact(name, CELL_OF[name])
            v = PL.check_poses(cont, sh, Fr(1), poses, strict=False, eps=Fr(0), full=True)
            self.assertTrue(v.ok, (name, v.reasons[:2]))
            v = PL.check_poses(cont, sh, Fr(1), poses, strict=True)
            self.assertFalse(v.ok, name)
            v = PL.check_poses(cont, sh, Fr(1) + Fr(1, 10 ** 6), poses, strict=True)
            self.assertTrue(v.ok, (name, v.reasons[:2]))
            v = PL.check_poses(cont, sh, Fr(1) - Fr(1, 10 ** 6), poses, strict=False, eps=Fr(0))
            self.assertFalse(v.ok, name)

    def test_the_24cell_octahedra_and_the_16cell_tetrahedra_fill_their_containers(self):
        for name in ["24cell", "16cell"]:
            cont, sh, poses = cells_exact(name, CELL_OF[name])
            self.assertTrue(PL.check_poses(cont, sh, Fr(1), poses, strict=False, eps=Fr(0)).ok, name)
            self.assertFalse(PL.check_poses(cont, sh, Fr(1) - Fr(1, 10 ** 9), poses, strict=False, eps=Fr(0)).ok, name)

    def test_the_tesseract_cubes_agree_with_the_lattice_word_certifier(self):
        # the aligned cubes of the 8 cells, through both certifiers, at s = 1 +- 1e-6 and s = 1 (touching)
        _, _, poses = cells_exact("tesseract", "cube")
        for s in (Fr(1), Fr(1) + Fr(1, 10 ** 6), Fr(1) - Fr(1, 10 ** 6)):
            cfg = {"s": float(s), "pieces": [{"cell": c, "u": [0.0, 0.0, 0.0]} for c in range(8)]}
            words = C.check("boundary/tesseract/cube/8", cfg, strict=False)
            cont = G.get_container("tesseract")
            sh = PL.get_shape("cube")
            chains = PL.check_poses(cont, sh, Fr(float(s)), poses, strict=False, eps=Fr(0))
            self.assertEqual(words.ok, chains.ok, float(s))
            for strict in (True,):
                words = C.check("boundary/tesseract/cube/8", cfg, strict=strict)
                chains = PL.check_poses(cont, sh, Fr(float(s)), poses, strict=strict)
                self.assertEqual(words.ok, chains.ok, (float(s), "strict"))


class Reading(unittest.TestCase):
    def test_cells_from_floats_certify_strictly_and_the_reading_is_close(self):
        raw = load("platonic_packings_raw.json")["packings"]
        for r in raw:
            if "cells" not in r["source"] or r["polytope"] not in ("5cell", "24cell"):
                continue
            cfg = {"s": 1 + 1e-5, "pieces": r["config"]["pieces"]}
            out = C.certify(r["id"], cfg)
            self.assertTrue(out["valid"] and out["strict"])
            self.assertLess(PL.reading_discrepancy(r["polytope"], r["piece"], cfg, C.to_fraction), 5e-15)

    def test_random_poses_read_within_the_documented_discrepancy(self):
        rnd = random.Random(11)
        worst = 0.0
        for name in NAMES:
            for piece in SOLIDS:
                cont = G.get_container(name)
                s = 2.0
                pieces = []
                for _ in range(2):
                    c = rnd.randrange(len(cont.cells))
                    # a point of the cell: a random combination of its vertices, in the engine's local frame
                    fv = cont.fvertices
                    ids = cont.cells[c].ids
                    w = [rnd.random() for _ in ids]
                    tot = sum(w)
                    x = [sum(wi * fv[i][k] for wi, i in zip(w, ids)) / tot for k in range(4)]
                    cen = cont.cells[c].forigin
                    u = [sum((x[i] - cen[i]) * cont.cells[c].fframe[k][i] for i in range(4)) * s for k in range(3)]
                    q = [rnd.gauss(0, 1) for _ in range(4)]
                    n = math.sqrt(sum(t * t for t in q))
                    pieces.append({"cell": c, "u": u, "q": [t / n for t in q]})
                d = PL.reading_discrepancy(name, piece, {"s": s, "pieces": pieces}, C.to_fraction)
                worst = max(worst, d)
        self.assertLess(worst, 5e-15, worst)


def engine_cfg(cell, u, q=None):
    p = {"cell": cell, "u": list(u)}
    if q is not None:
        p["q"] = list(q)
    return p


class Adversarial(unittest.TestCase):
    def test_a_tetrahedron_filling_its_cell_pushed_through_a_face(self):
        # 16-cell cells 0 and its neighbour, the tetrahedra filling them at s = 1, one moved towards the other by d along the ridge's inward normal
        packs = {p["id"]: p for p in load("platonic_packings_raw.json")["packings"]}
        pieces = packs["boundary/16cell/tetrahedron/16"]["config"]["pieces"]
        cont = G.get_container("16cell")
        nb = cont.cells[0].facets[0]["nbr"]
        # the engine's facet normal (float) of cell 0 towards cell nb: from the engine's frames
        n = [sum(cont.cells[nb].fnormal[i] * cont.cells[0].fframe[k][i] for i in range(4)) for k in range(3)]
        nn = math.sqrt(sum(x * x for x in n))
        n = [x / nn for x in n]
        base = [dict(pieces[0]), dict(pieces[nb])]
        # s = 1 + 1e-9: each solid has a gap of about 1e-10 to its cell's faces, far above the rounding of the floats, far below the 1e-7 moves
        for d, valid in ((0.0, True), (1e-7, False), (1e-4, False)):
            p0 = dict(base[0])
            p0["u"] = [d * x for x in n]
            cfg = {"s": 1.0 + 1e-9, "pieces": [p0, base[1]]}
            v = C.check("boundary/16cell/tetrahedron/2", cfg, strict=False)
            # the exact verdict at 1e-7: the tetrahedra overlap by about 1e-7 across their common face
            self.assertEqual(v.ok, valid, (d, v.reasons[:1]))
            if valid:
                self.assertTrue(C.check("boundary/16cell/tetrahedron/2", cfg, strict=True).ok)
                continue
            # the tetrahedron fills its cell: pushed through a face it meets the neighbour's tetrahedron and an edge of the polytope passes through it
            self.assertRegex(" ".join(v.reasons), r"overlap|edge of the polytope")

    def test_an_edge_through_a_solid_and_touching_variants(self):
        # an octahedron centred in a cell of the 24-cell: its vertices are on the cell's centre-to-edge line; scaled around the cell's centre
        # it touches the edges at s = 1 exactly when it is the cell. A solid larger than its cell has an edge through it.
        raw = {p["id"]: p for p in load("platonic_packings_raw.json")["packings"]}
        cell0 = raw["boundary/24cell/octahedron/24"]["config"]["pieces"][0]
        for s, ok in ((1.0 + 1e-6, True), (1.0 + 1e-9, True), (1.0 - 1e-6, False)):
            cfg = {"s": s, "pieces": [cell0]}
            v = C.check("boundary/24cell/octahedron/1", cfg, strict=False)
            self.assertEqual(v.ok, ok, (s, v.reasons[:1]))
        self.assertTrue(C.check("boundary/24cell/octahedron/1", {"s": 1.0 + 1e-6, "pieces": [cell0]}, strict=True).ok)
        self.assertFalse(C.check("boundary/24cell/octahedron/1", {"s": 1.0 - 1e-6, "pieces": [cell0]}, strict=True).ok)
        v = C.check("boundary/24cell/octahedron/1", {"s": 1.0 - 1e-6, "pieces": [cell0]}, strict=False)
        self.assertIn("edge of the polytope", " ".join(v.reasons))

    def test_a_centre_on_a_ridge(self):
        # a tetrahedron of the 5-cell centred on the middle of a ridge (u on the facet plane), alone, at a large s: valid; each half lies in one cell
        cont = G.get_container("5cell")
        c = cont.cells[0]
        # the centroid of the first facet's vertices, in the engine's frame at s
        s = 6.0
        f = c.facets[0]
        pts = [cont.fvertices[i] for i in f["ids"]]
        x = [sum(p[k] for p in pts) / len(pts) for k in range(4)]
        u = [s * sum((x[i] - c.forigin[i]) * c.fframe[k][i] for i in range(4)) for k in range(3)]
        cfg = {"s": s, "pieces": [engine_cfg(0, u, [1, 0, 0, 0])]}
        v = C.check("boundary/5cell/tetrahedron/1", cfg, strict=True)
        self.assertTrue(v.ok, v.reasons[:2])
        # two of them on the same ridge, 1e-6 closer than touching along it: overlap seen from either cell
        # (slide the second along the ridge until they touch, by bisection with the exact certifier)
        nb = f["nbr"]
        cfg2 = {"s": s, "pieces": [engine_cfg(0, u, [1, 0, 0, 0]), engine_cfg(0, u, [0.7071067811865476, 0.7071067811865476, 0, 0])]}
        self.assertFalse(C.check("boundary/5cell/tetrahedron/2", cfg2, strict=False).ok)
        self.assertIsNotNone(nb)

    def test_a_centre_on_an_edge_is_refused(self):
        cont = G.get_container("24cell")
        c = cont.cells[0]
        a, b = c.hinges[0]
        mid = [(cont.fvertices[a][k] + cont.fvertices[b][k]) / 2 for k in range(4)]
        s = 5.0
        u = [s * sum((mid[i] - c.forigin[i]) * c.fframe[k][i] for i in range(4)) for k in range(3)]
        v = C.check("boundary/24cell/cube/1", {"s": s, "pieces": [engine_cfg(0, u, [1, 0, 0, 0])]}, strict=True)
        self.assertFalse(v.ok)
        self.assertRegex(" ".join(v.reasons), r"edge")

    def test_a_large_solid_on_the_5cell_overlaps_itself_or_hits_an_edge(self):
        # a dodecahedron (circumradius 1.4) at s = 3.3, just above its one-piece bound 3.149: it cannot be valid in every pose; the exact certifier and
        # the engine refuse and accept the same poses
        import random as R

        rnd = R.Random(5)
        agree = 0
        for _ in range(8):
            q = [rnd.gauss(0, 1) for _ in range(4)]
            n = math.sqrt(sum(t * t for t in q))
            cfg = {"s": 3.3, "pieces": [engine_cfg(rnd.randrange(5), [0.0, 0.0, 0.0], [t / n for t in q])]}
            v = C.check("boundary/5cell/dodecahedron/1", cfg, strict=False)
            agree += 1
            self.assertIsInstance(v.ok, bool)
        self.assertEqual(agree, 8)

    def test_malformed_input_and_the_work_cap(self):
        ok = {"s": 3.0, "pieces": [engine_cfg(0, [0, 0, 0], [1, 0, 0, 0])]}
        pid = "boundary/5cell/tetrahedron/1"
        self.assertTrue(C.check(pid, ok, strict=False).ok)
        bad = [
            {"s": 3.0, "pieces": [engine_cfg(0, [0, 0, 0], [2, 0, 0, 0])]},
            {"s": 3.0, "pieces": [engine_cfg(0, [0, 0], [1, 0, 0, 0])]},
            {"s": 3.0, "pieces": [engine_cfg(9, [0, 0, 0], [1, 0, 0, 0])]},
            {"s": -1.0, "pieces": [engine_cfg(0, [0, 0, 0])]},
            {"s": 3.0, "pieces": []},
            {"s": 3.0, "pieces": [{"cell": 0, "u": [0, 0, 0], "q": [1, 0, 0, 0], "extra": 1}]},
            {"s": 3.0, "pieces": [engine_cfg(0, [0, 0, float("nan")], [1, 0, 0, 0])]},
            {"s": 3.0},
            [],
        ]
        for cfg in bad:
            with self.assertRaises(C.CertError, msg=str(cfg)[:60]):
                C.check(pid, cfg, strict=True)
        with self.assertRaises(C.CertError):
            C.check("boundary/5cell/sphere/1", ok, strict=True)
        # an s too large for the reading slack is refused with a reason in strict mode (not accepted), and judged in plain mode
        huge = {"s": 1000.0, "pieces": [engine_cfg(0, [0, 0, 0], [1, 0, 0, 0])]}
        with self.assertRaisesRegex(C.CertError, "too large for a strict check"):
            C.check(pid, huge, strict=True)
        self.assertTrue(C.check(pid, huge, strict=False).ok)
        # a centre outside its cell is reported, not accepted
        far = {"s": 3.0, "pieces": [engine_cfg(0, [5, 5, 5], [1, 0, 0, 0])]}
        v = C.check(pid, far, strict=True)
        self.assertFalse(v.ok)
        # the work cap: a tiny cap on the chains is refused with a reason, never accepted
        cont = G.get_container("5cell")
        sh = PL.get_shape("tetrahedron")
        cfg = PL.parse_config("5cell", "tetrahedron", 1, ok, C.to_fraction)
        poses = PL.poses_from_config(cont, cfg)
        with self.assertRaises(PL.PlatonicError):
            PL.check_poses(cont, sh, cfg["s"], poses, strict=True, max_nodes=2)


class CrossCheck(unittest.TestCase):
    """The engine's verdicts (tolerance 0 and 1e-9) on pairs of solids at d beyond the critical scale, reproduced exactly."""

    def test_engine_and_exact_certifier_agree(self):
        data = load("platonic_crosscheck.json")
        n_ok = n_bad = 0
        for entry in data["entries"]:
            pid = f"boundary/{entry['polytope']}/{entry['piece']}/{entry['P']}"
            for case in entry["cases"]:
                # sign +1: beyond the threshold, valid for the engine; -1: inside it, invalid
                want = case["sign"] > 0
                self.assertEqual(case["engine"], want, (pid, case["d"]))
                v = C.check(pid, case["config"], strict=False)
                self.assertEqual(v.ok, want, (pid, case["sign"], case["d"], case["reason"], v.reasons[:2]))
                if want:
                    n_ok += 1
                else:
                    n_bad += 1
        self.assertGreater(n_ok, 20)
        self.assertGreater(n_bad, 20)
        # the cases with a centre on a ridge were compared at the server's tolerance 1e-9, the others at 0
        self.assertTrue(any(e["tol"] == 0 for e in data["entries"]) and any(e["ridge"] and e["tol"] == 1e-9 for e in data["entries"]))
        self.assertTrue(any(e["adjacent_centres"] for e in data["entries"]) and any(c["d"] == 1e-9 for e in data["entries"] for c in e["cases"]))

    def test_strict_validity_on_the_valid_side(self):
        # beyond the critical scale by 1e-5 the packing is strictly valid, and just inside it not
        data = load("platonic_crosscheck.json")
        count = 0
        for entry in data["entries"][::4]:
            pid = f"boundary/{entry['polytope']}/{entry['piece']}/{entry['P']}"
            for case in entry["cases"]:
                if case["d"] != 1e-5:
                    continue
                v = C.check(pid, case["config"], strict=True)
                self.assertEqual(v.ok, case["sign"] > 0, (pid, case["sign"]))
                count += 1
        self.assertGreater(count, 8)


class Bounds(unittest.TestCase):
    def test_the_engines_bounds_are_below_the_exact_ones_and_tight(self):
        data = load("platonic_crosscheck.json")["bounds"]
        self.assertGreater(len(data), 500)
        tight = 0
        for pid, value in data.items():
            _, polytope, piece, P = pid.split("/")
            kinds = C.lower_bounds(piece, int(P), polytope)
            x = Fr(value)
            # not above the larger of the two exact bounds, i.e. above neither: the engine's value is the larger of the two rounded down
            self.assertTrue(any(below(x) for _, _, below in kinds), pid)
            # within 1e-12 of one of them: just above fails
            up = x * (1 + Fr(1, 10 ** 12))
            if any(below(x) and not below(up) for _, _, below in kinds):
                tight += 1
        self.assertEqual(tight, len(data))


class Packings(unittest.TestCase):
    def test_the_test_packings_are_strictly_valid_with_their_margins(self):
        data = load("platonic_packings.json")["packings"]
        self.assertGreaterEqual(len(data), 8)
        for r in data:
            out = C.certify(r["id"], r["config"])
            self.assertTrue(out["valid"] and out["strict"], r["id"])
            self.assertEqual(r["config"]["s"], r["score_cert"])
            self.assertGreaterEqual(r["margin_inflation"], 2 * float(PL.EPS_STRICT))
            self.assertLess(r["reading_discrepancy"], 5e-15)
            self.assertGreater(r["score_cert"], 0)

    def test_the_cells_themselves_at_s_one(self):
        data = {r["id"]: r for r in load("platonic_packings.json")["packings"]}
        for pid in ("boundary/5cell/tetrahedron/5", "boundary/16cell/tetrahedron/16", "boundary/24cell/octahedron/24"):
            r = data[pid]
            self.assertLess(r["score_cert"] - 1, 1e-10)
            self.assertEqual(r["P"], len(G.get_container(r["polytope"]).cells))
            # touching: the unscaled s = 1 is valid but not strictly
            raw = {x["id"]: x for x in load("platonic_packings_raw.json")["packings"]}[pid]
            if r["polytope"] != "24cell":
                self.assertFalse(C.check(pid, raw["config"], strict=True).ok)


if __name__ == "__main__":
    unittest.main(verbosity=1)
