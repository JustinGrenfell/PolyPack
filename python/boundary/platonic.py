"""Exact certifier for Platonic pieces on the boundary of a regular 4-polytope ("3D on 4D"; BOUNDARY.md section 4C).

Problem boundary/<polytope>/<piece>/<P>, piece one of tetrahedron, cube, octahedron, dodecahedron, icosahedron, on the tesseract (every
piece but its axis-aligned cube, which certify.py's lattice words judge), 5-cell, 16-cell, 24-cell, 120-cell or 600-cell: P unit-edge
solids, freely rotated, on the boundary of the polytope of edge s; the score is s. The conventions are those of
packages/boundary/src/types.ts and general.py: a configuration is {s, pieces: [{cell, u, q?}]}, u the piece's centre in the local frame of
its anchor cell at the configuration's s, q the unit quaternion (w, x, y, z) that rotates the standard solid (packages/boundary/src/pieces.ts)
in that frame.

What is decided, exactly (no floating point decides anything; floats only choose which exact test to try first and prune with rigorous
error bounds):

  * Every number is an element of Q(sqrt 2, sqrt 5) (qmulti.QM), or a rational. The container is general.Container (its charts, gluing
    rotations and chain enumerator are reused unchanged, in Q(sqrt 5)); the pieces' coordinates add sqrt 2 and sqrt 5 as needed.
  * The piece is placed through an exact orthonormal frame of its cell's hyperplane. Cell by cell (cell_frame) three mutually orthogonal
    vectors f_k of equal squared length m are built from the cell's own vertices (the tetrahedron's circumscribed-cube axes, the octahedron's
    diameters, the cube's edges, the dodecahedron's inscribed-cube edges), so that e0 / sqrt(m) = kappa is 1/sqrt 2 (tetrahedra, octahedra),
    1/phi (dodecahedra) or 1 (cubes); the base polytope has edge e0 and the piece, of unit edge in a polytope of edge s, has edge e0 / s there.
    The piece's vertex v is  x + (kappa / s) (R v_std)_k f_k, k = 1..3, x the centre, R an exact rational rotation (a rational quaternion
    gives one), v_std the standard solid's exact vertex.
  * Reading a configuration. The engine's floats mean the ambient point s * origin + sum u_k frame_k and the rotation of the standard solid
    by q in the engine's float frame. The centre is read as in general.py (an exact point of the cell's hyperplane, dyadic weights of four vertices).
    The rotation: the engine's frame and the exact frame span the same 3-space, so the engine's rotation in the exact frame is O R_e
    (O = F_x^T F_e, both frames in the engine's float coordinates); that float matrix is turned into a quaternion, rounded to a dyadic grid of
    2^-60 and made an exact rational rotation. The discrepancy of the whole reading (centres and rotations) is measured by reading_discrepancy().
  * Strict validity is decided on the piece inflated about its centre by the factor 1 + eps (EPS_STRICT): the inflated pieces must be valid
    (interiors disjoint, no edge of the polytope through an interior, embedded). Then the closed pieces are strictly valid with room: each is at
    distance at least eps * inradius from the others and from the 1-skeleton, which is more than the reading discrepancy (documented, 4C.3), so
    "strictly valid" binds the configuration as the engine reads it. eps = 0 decides plain validity of the exact poses.
  * Development. For each piece and each closed cell that contains its centre, the chains of cells within the piece's circumradius are
    enumerated (general.enumerate_chains, exact). Along a chain C ending in cell c the part of the piece in c is X ∩ c ∩ W(C), mapped into
    c's own chart. It is computed by exact clipping of the piece's vertex set by the cell's facet planes and the window's planes, tracking which
    planes each vertex lies on (so the part's edges are known without ever ordering anything). Parts without interior are dropped.
  * The decisions. (1) No edge of the polytope meets the interior of a developed piece: every hinge of every node's cell, clipped to the closed
    window, against the open interior of the piece (an interval test per face, exact). (2) In every cell the parts are pairwise interior-disjoint,
    parts of one piece included: the exact separating-axis test on face normals and edge-edge cross products of the two parts; a pair is accepted
    only with an exact separating plane, and a pair without one is reported as overlapping.

Usage: imported by certify.py. python3 python/boundary/platonic.py boundary/<polytope>/<piece>/<P> <config.json> certifies one configuration.
"""

import itertools
import math
from fractions import Fraction

import general as G
from general import Q5, cross, inv3, matvec, vadd, vdot, vscale, vsub
from qmulti import QM

PIECES = ("tetrahedron", "cube", "octahedron", "dodecahedron", "icosahedron")

# Strict validity is decided on pieces inflated by this factor (module docstring). 2^-40 is 9.1e-13: the engine's reading is within a
# few 1e-15 s of the exact poses, and the fixtures' own scale (1 + k 1e-12) leaves contact gaps of about 1e-12 times the piece's size.
EPS_STRICT = Fraction(1, 2 ** 40)

MAX_NODES = 40000  # chains per piece and run before the check is refused
MAX_S = Fraction(10) ** 4
QUAT_BITS = 60


class PlatonicError(Exception):
    """Malformed input or a configuration beyond the certifier's work cap (certify.CertError at the front end)."""


# ------------------------------------------------------------------ small helpers


def qs(x):
    """-1, 0 or 1: the exact sign of a QM, Q5, Fraction or int."""
    if isinstance(x, (QM, Q5)):
        return x.sign()
    return (x > 0) - (x < 0)


def qm(x):
    return QM.of(x)


def qvec(v):
    return [QM.of(x) for x in v]


def popcount(x):
    return bin(x).count("1")


def fl(x):
    """A binary64 value of x (a QM, Q5, Fraction, int)."""
    if isinstance(x, QM):
        return x.fl()[0]
    return float(x)


# ------------------------------------------------------------------ the standard solids (exact)

_SHAPES = {}


def _signs(base):
    out = []
    for m in range(8):
        v = tuple(-x if (m >> k) & 1 else x for k, x in enumerate(base))
        if v not in out:
            out.append(v)
    return out


def _std_vertices(name):
    half = Fraction(1, 2)
    phi = QM(half, 0, half)  # (1 + sqrt 5) / 2
    zero = QM(0)
    if name == "tetrahedron":
        k = QM(0, Fraction(1, 4))  # 1 / (2 sqrt 2)
        return [[k * x for x in p] for p in ((1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1))]
    if name == "cube":
        return [[QM(half if (v >> k) & 1 else -half) for k in range(3)] for v in range(8)]
    if name == "octahedron":
        r = QM(0, half)  # 1 / sqrt 2
        out = []
        for k in range(3):
            for sg in (1, -1):
                p = [zero, zero, zero]
                p[k] = r * sg
                out.append(p)
        return out
    if name == "dodecahedron":
        a = phi * half
        b = phi * phi * half
        out = []
        for sx in (1, -1):
            for sy in (1, -1):
                for sz in (1, -1):
                    out.append([a * sx, a * sy, a * sz])
        seen = []
        for base in ((zero, QM(half), b), (QM(half), b, zero), (b, zero, QM(half))):
            for sg in itertools.product((1, -1), repeat=3):
                p = [base[0] * sg[0], base[1] * sg[1], base[2] * sg[2]]
                key = tuple(p)
                if key not in seen:
                    seen.append(key)
                    out.append(p)
        return out
    if name == "icosahedron":
        a = phi * half
        out = []
        seen = []
        for base in ((zero, QM(half), a), (QM(half), a, zero), (a, zero, QM(half))):
            for sg in itertools.product((1, -1), repeat=3):
                p = [base[0] * sg[0], base[1] * sg[1], base[2] * sg[2]]
                key = tuple(p)
                if key not in seen:
                    seen.append(key)
                    out.append(p)
        return out
    raise PlatonicError(f"unknown piece {name!r}")


class Shape:
    """A standard solid of unit edge centred at the origin: exact vertices, faces (vertex sets with a generating triple), the exact squared
    circumradius and inradius."""


def get_shape(name):
    if name in _SHAPES:
        return _SHAPES[name]
    V = _std_vertices(name)
    n = len(V)
    sh = Shape()
    sh.name, sh.V = name, V
    # unit edge: the least squared distance between vertices is exactly 1, and the number of such pairs is the literature's
    d2 = [[vdot(vsub(V[i], V[j]), vsub(V[i], V[j])) for j in range(n)] for i in range(n)]
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n) if d2[i][j] == 1]
    counts = {"tetrahedron": 6, "cube": 12, "octahedron": 12, "dodecahedron": 30, "icosahedron": 30}
    if len(pairs) != counts[name]:
        raise PlatonicError(f"{name}: edge count {len(pairs)} differs from the literature")
    sh.edges = pairs
    fV = [[fl(c) for c in v] for v in V]
    faces = {}
    for a, b, c in itertools.combinations(range(n), 3):
        e1 = [fV[b][k] - fV[a][k] for k in range(3)]
        e2 = [fV[c][k] - fV[a][k] for k in range(3)]
        N = [e1[1] * e2[2] - e1[2] * e2[1], e1[2] * e2[0] - e1[0] * e2[2], e1[0] * e2[1] - e1[1] * e2[0]]
        if math.sqrt(sum(x * x for x in N)) < 1e-9:
            continue
        side = [sum(N[k] * (fV[i][k] - fV[a][k]) for k in range(3)) for i in range(n)]
        if any(x > 1e-9 for x in side) and any(x < -1e-9 for x in side):
            continue
        Ne = cross(vsub(V[b], V[a]), vsub(V[c], V[a]))
        sg = [qs(vdot(Ne, vsub(V[i], V[a]))) for i in range(n)]
        if 1 in sg and -1 in sg:
            continue
        on = tuple(i for i in range(n) if sg[i] == 0)
        if on not in faces:
            faces[on] = (a, b, c)
    sh.faces = [{"verts": on, "tri": tri} for on, tri in sorted(faces.items())]
    fcount = {"tetrahedron": 4, "cube": 6, "octahedron": 8, "dodecahedron": 12, "icosahedron": 20}
    if len(sh.faces) != fcount[name]:
        raise PlatonicError(f"{name}: face count {len(sh.faces)} differs from the literature")
    sh.vmask = [sum(1 << k for k, f in enumerate(sh.faces) if i in f["verts"]) for i in range(n)]
    sh.circ2 = _as_q5(vdot(V[0], V[0]))
    # inradius^2: the squared distance from the origin to a face plane, N.p / |N| with p a vertex of the face
    a, b, c = sh.faces[0]["tri"]
    Ne = cross(vsub(V[b], V[a]), vsub(V[c], V[a]))
    h = vdot(Ne, V[a])
    sh.in2 = _as_q5(h * h / vdot(Ne, Ne))
    sh.volume2 = None
    _SHAPES[name] = sh
    return sh


def _as_q5(x):
    """A QM known to lie in Q(sqrt 5) as a Q5."""
    x = QM.of(x)
    if x.n1 != 0 or x.n3 != 0:
        raise PlatonicError("expected a number of Q(sqrt 5)")
    return Q5(Fraction(x.n0, x.d), Fraction(x.n2, x.d))


# ------------------------------------------------------------------ frames of the cells (exact)


class CellFrame:
    pass


def cell_frame(cont, ci):
    """Three mutually orthogonal vectors of the hyperplane of cell ci, equal squared length m, from the cell's own vertices; kappa (a QM) with
    kappa^2 m = e0^2; phi_k = their chart images (Q5). Cached on the cell."""
    c = cont.cells[ci]
    if hasattr(c, "frame"):
        return c.frame
    V = cont.V
    P = [V[i] for i in c.ids]
    name = cont.name
    zero = cont.zero
    if len(P) == 4:  # tetrahedral cell: the axes of the circumscribed cube
        f = [
            vsub(vadd(P[0], P[1]), vadd(P[2], P[3])),
            vsub(vadd(P[0], P[2]), vadd(P[1], P[3])),
            vsub(vadd(P[0], P[3]), vadd(P[1], P[2])),
        ]
    elif len(P) == 6:  # octahedral cell: three diameters
        f = []
        used = set()
        for i in range(6):
            if i in used:
                continue
            for j in range(i + 1, 6):
                d = vsub(P[i], P[j])
                if vdot(d, d) == 2 * cont.e2 and j not in used:
                    f.append(d)
                    used |= {i, j}
                    break
        if len(f) != 3:
            raise PlatonicError("octahedral cell without three diameters")
    elif len(P) == 8:  # cubic cell: three edges from a vertex
        O = P[0]
        f = [vsub(p, O) for p in P[1:] if vdot(vsub(p, O), vsub(p, O)) == cont.e2]
        if len(f) != 3:
            raise PlatonicError("cubic cell without three edges at a vertex")
    elif len(P) == 20:  # dodecahedral cell: three edges of an inscribed cube (face diagonals, length phi e)
        phi = Q5(Fraction(1, 2), Fraction(1, 2))
        D2 = phi * phi * cont.e2
        f = None
        for i in range(20):
            cand = [vsub(P[j], P[i]) for j in range(20) if j != i and vdot(vsub(P[j], P[i]), vsub(P[j], P[i])) == D2]
            for tri in itertools.combinations(cand, 3):
                if all(vdot(tri[a], tri[b]) == 0 for a in range(3) for b in range(a + 1, 3)):
                    f = list(tri)
                    break
            if f:
                break
        if not f:
            raise PlatonicError("dodecahedral cell without an inscribed-cube triad")
    else:
        raise PlatonicError("unknown cell")
    m = vdot(f[0], f[0])
    for a in range(3):
        if vdot(f[a], f[a]) != m:
            raise PlatonicError("frame vectors of unequal length")
        for b in range(a + 1, 3):
            if vdot(f[a], f[b]) != 0:
                raise PlatonicError("frame vectors not orthogonal")
    # kappa = e0 / sqrt m: the cell kind fixes it (module docstring); the exact check kappa^2 m = e0^2 guards it
    if len(P) in (4, 6):
        kappa = QM(0, Fraction(1, 2))
    elif len(P) == 8:
        kappa = QM(1)
    else:
        kappa = QM(Fraction(-1, 2), 0, Fraction(1, 2))  # 1 / phi = (sqrt 5 - 1) / 2
    if kappa * kappa * QM.of(m) != QM.of(cont.e2):
        raise PlatonicError("kappa^2 m differs from e0^2")
    fr = CellFrame()
    fr.f, fr.m, fr.kappa = f, m, kappa
    fr.phi = [cont.chart_lin(c, fk) for fk in f]  # chart image of each frame vector (Q5 3-vectors)
    c.frame = fr
    return fr


# ------------------------------------------------------------------ the engine's frames in floats


def _matmul3(a, b):
    return [[sum(a[i][k] * b[k][j] for k in range(3)) for j in range(3)] for i in range(3)]


def _det3f(m):
    return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])


def quat_matrix_float(q):
    """The twin of packages/boundary/src/linalg.ts quaternionToMatrix: the rotation matrix of the normalised quaternion (w, x, y, z)."""
    n = math.hypot(q[0], q[1], q[2], q[3])
    w, x, y, z = (q[0] / n, q[1] / n, q[2] / n, q[3] / n)
    return [
        [1 - 2 * (y * y + z * z), 2 * (x * y - w * z), 2 * (x * z + w * y)],
        [2 * (x * y + w * z), 1 - 2 * (x * x + z * z), 2 * (y * z - w * x)],
        [2 * (x * z - w * y), 2 * (y * z + w * x), 1 - 2 * (x * x + y * y)],
    ]


def matrix_quaternion(m):
    """Unit quaternion (w, x, y, z) of a proper rotation matrix (list of rows)."""
    tr = m[0][0] + m[1][1] + m[2][2]
    if tr > 0:
        S = math.sqrt(tr + 1) * 2
        w, x, y, z = S / 4, (m[2][1] - m[1][2]) / S, (m[0][2] - m[2][0]) / S, (m[1][0] - m[0][1]) / S
    elif m[0][0] > m[1][1] and m[0][0] > m[2][2]:
        S = math.sqrt(1 + m[0][0] - m[1][1] - m[2][2]) * 2
        w, x, y, z = (m[2][1] - m[1][2]) / S, S / 4, (m[0][1] + m[1][0]) / S, (m[0][2] + m[2][0]) / S
    elif m[1][1] > m[2][2]:
        S = math.sqrt(1 + m[1][1] - m[0][0] - m[2][2]) * 2
        w, x, y, z = (m[0][2] - m[2][0]) / S, (m[0][1] + m[1][0]) / S, S / 4, (m[1][2] + m[2][1]) / S
    else:
        S = math.sqrt(1 + m[2][2] - m[0][0] - m[1][1]) * 2
        w, x, y, z = (m[1][0] - m[0][1]) / S, (m[0][2] + m[2][0]) / S, (m[1][2] + m[2][1]) / S, S / 4
    n = math.sqrt(w * w + x * x + y * y + z * z)
    return [w / n, x / n, y / n, z / n]


def rational_rotation(q):
    """The exact rotation matrix (3 x 3 Fractions) of the rational quaternion obtained by rounding q to the grid 2^-QUAT_BITS."""
    iq = [round(x * (1 << QUAT_BITS)) for x in q]
    w, x, y, z = iq
    N = w * w + x * x + y * y + z * z
    if N == 0:
        raise PlatonicError("zero quaternion")
    rows = [
        [w * w + x * x - y * y - z * z, 2 * (x * y - w * z), 2 * (x * z + w * y)],
        [2 * (x * y + w * z), w * w - x * x + y * y - z * z, 2 * (y * z - w * x)],
        [2 * (x * z - w * y), 2 * (y * z + w * x), w * w - x * x - y * y + z * z],
    ]
    return [[Fraction(v, N) for v in r] for r in rows]


def engine_frame_images(cont, ci):
    """The cell's exact frame vectors as unit vectors of the engine's float ambient coordinates: F_x[k] (4 floats each), by the affine map that
    takes the cell's exact vertices to the engine's unit-edge float vertices."""
    c = cont.cells[ci]
    fr = cell_frame(cont, ci)
    fv = cont.fvertices
    e0 = math.sqrt(float(cont.e2))
    scale = e0 / math.sqrt(float(fr.m))
    rows = []
    for k in range(3):
        lam = [float(x) for x in fr.phi[k]]
        v = [0.0] * 4
        for j in range(3):
            d = [fv[c.tri[j]][i] - fv[c.O_idx][i] for i in range(4)]
            for i in range(4):
                v[i] += lam[j] * d[i]
        rows.append([x * scale for x in v])
    # orthonormality (a guard on the whole chain of conventions)
    for a in range(3):
        for b in range(3):
            g = sum(rows[a][i] * rows[b][i] for i in range(4))
            if abs(g - (1.0 if a == b else 0.0)) > 1e-9:
                raise PlatonicError("the exact frame is not orthonormal in the engine's coordinates")
    return rows


# ------------------------------------------------------------------ reading a configuration


class PieceReading:
    """One piece, read exactly: its anchor, centre (chart coordinates in every closed cell containing it), the exact rotation in the cell's
    frame and the sign of the third frame vector."""


def parse_config(polytope, piece, P, config, to_fraction):
    cont = G.get_container(polytope)
    if piece not in PIECES:
        raise PlatonicError(f"unknown piece {piece!r}")
    if not isinstance(config, dict):
        raise PlatonicError("the configuration must be an object with s and pieces")
    s = to_fraction(config.get("s"), "s")
    if not 0 < s <= MAX_S:
        raise PlatonicError("s must be greater than 0 and at most 1e4")
    pieces = config.get("pieces")
    if not isinstance(pieces, list):
        raise PlatonicError("the configuration needs a list of pieces")
    if len(pieces) != P:
        raise PlatonicError(f"the configuration has {len(pieces)} pieces; the problem needs {P}")
    out = []
    C = len(cont.cells)
    for n, p in enumerate(pieces):
        if not isinstance(p, dict):
            raise PlatonicError(f"piece {n} must be an object with cell, u and q")
        extra = set(p) - {"cell", "u", "q"}
        if extra:
            raise PlatonicError(f"piece {n} has unknown fields {sorted(extra)}")
        c = p.get("cell")
        if isinstance(c, bool) or not isinstance(c, int) or not 0 <= c < C:
            raise PlatonicError(f"piece {n}: cell must be a whole number from 0 to {C - 1}")
        u = p.get("u")
        if not isinstance(u, list) or len(u) != 3:
            raise PlatonicError(f"piece {n}: u must be 3 numbers")
        uu = [to_fraction(x, f"piece {n} u[{k}]") for k, x in enumerate(u)]
        if any(abs(x) > MAX_S for x in uu):
            raise PlatonicError(f"piece {n}: u is out of range")
        q = p.get("q")
        if q is None:
            qq = [Fraction(1), Fraction(0), Fraction(0), Fraction(0)]
        else:
            if not isinstance(q, list) or len(q) != 4:
                raise PlatonicError(f"piece {n}: q must be 4 numbers (a unit quaternion w, x, y, z)")
            qq = [to_fraction(x, f"piece {n} q[{k}]") for k, x in enumerate(q)]
            nrm = math.sqrt(float(sum(x * x for x in qq)))
            if not abs(nrm - 1) <= 1e-9:
                raise PlatonicError(f"piece {n}: q must be a unit quaternion (norm within 1e-9 of 1)")
        out.append({"cell": c, "u": uu, "q": qq})
    return {"s": s, "pieces": out}


def read_rotation(cont, ci, q):
    """(R, sign3): the exact rational rotation of the piece in the exact frame of cell ci that the engine's float quaternion q means, and the sign
    (+1 or -1) of the third frame vector that makes the two frames' handedness agree."""
    c = cont.cells[ci]
    Fx = engine_frame_images(cont, ci)
    Fe = c.fframe
    O = [[sum(Fx[a][i] * Fe[b][i] for i in range(4)) for b in range(3)] for a in range(3)]
    sign3 = 1
    if _det3f(O) < 0:
        sign3 = -1
        O[2] = [-x for x in O[2]]
    Re = quat_matrix_float([float(x) for x in q])
    Rx = _matmul3(O, Re)
    return rational_rotation(matrix_quaternion(Rx)), sign3


class Pose:
    """An exact piece pose in one anchor: cell, centre chart, rotation R (3x3 Fractions or QM, in the frame of the cell with the third vector
    multiplied by sign3)."""

    def __init__(self, cell, lam, R, sign3=1, point=None):
        self.cell, self.lam, self.R, self.sign3, self.point = cell, lam, R, sign3, point


def poses_from_config(cont, cfg, points=None):
    """The exact poses of a parsed configuration: centres as general.Problem reads them (an exact point of the cell's hyperplane), rotations by
    read_rotation."""
    prob = G.Problem(cont, {"s": cfg["s"], "pieces": [{"cell": p["cell"], "u": p["u"]} for p in cfg["pieces"]]}, points)
    poses = []
    for b, p in zip(prob.balls, cfg["pieces"]):
        R, sign3 = read_rotation(cont, b.cell, p["q"])
        poses.append(Pose(b.cell, b.cells[b.cell], R, sign3, point=b.x))
        poses[-1].closed = dict(b.cells)
    return poses


# ------------------------------------------------------------------ pieces in charts


def piece_vertices_in_anchor(cont, shape, s, pose, eps):
    """The inflated piece's vertices in the chart of the pose's anchor cell (QM 3-vectors), and its centre."""
    fr = cell_frame(cont, pose.cell)
    kappa = fr.kappa
    scale = kappa * QM.of(Fraction(1) / s * (1 + eps))
    Phi = fr.phi
    sg = pose.sign3
    lam0 = qvec(pose.lam)
    out = []
    for v in shape.V:
        rv = [vdot([QM.of(x) for x in pose.R[i]], v) * scale for i in range(3)]
        if sg < 0:
            rv[2] = -rv[2]
        P = []
        for i in range(3):
            acc = lam0[i]
            for k in range(3):
                acc = acc + rv[k] * QM.of(Phi[k][i])
            P.append(acc)
        out.append(P)
    return lam0, out


# ------------------------------------------------------------------ convex parts by incremental clipping


class Part:
    """A convex polytope in a cell's chart: exact vertices, for each vertex the bitmask of the planes it lies on, and the planes (N, h) with the
    polytope on the side N . y <= h."""

    __slots__ = ("V", "M", "planes", "fV", "box", "edges", "meta")


def affine_rank3(V):
    if len(V) < 4:
        return False
    o = V[0]
    d = [vsub(v, o) for v in V[1:]]
    for i in range(len(d)):
        if all(qs(x) == 0 for x in d[i]):
            continue
        for j in range(i + 1, len(d)):
            c = cross(d[i], d[j])
            if all(qs(x) == 0 for x in c):
                continue
            for k in range(j + 1, len(d)):
                if qs(vdot(c, d[k])) != 0:
                    return True
    return False


def adjacency(M):
    """Edges of the polytope from the vertex-plane incidence: vertices i, j are adjacent when at least two planes hold both and no other
    vertex lies on all the planes that hold both."""
    n = len(M)
    out = []
    for i in range(n):
        for j in range(i + 1, n):
            common = M[i] & M[j]
            if common & (common - 1) == 0:  # fewer than two common planes
                continue
            if sum(1 for k in range(n) if M[k] & common == common) == 2:
                out.append((i, j))
    return out


def clip(V, M, planes, N, h):
    """Clips the polytope (V, M, planes) to N . y <= h, exactly. Appends the plane to `planes` and returns the new (V, M), or None when no
    interior can remain (all vertices on or beyond the plane)."""
    d = [vdot(N, v) - h for v in V]
    sg = [qs(x) for x in d]
    if all(x <= 0 for x in sg):
        bit = 1 << len(planes)
        planes.append((N, h))
        return V, [m | (bit if sg[i] == 0 else 0) for i, m in enumerate(M)]
    if all(x >= 0 for x in sg):
        return None
    bit = 1 << len(planes)
    planes.append((N, h))
    newV, newM = [], []
    for i, v in enumerate(V):
        if sg[i] <= 0:
            newV.append(v)
            newM.append(M[i] | (bit if sg[i] == 0 else 0))
    for i, j in adjacency(M):
        if sg[i] * sg[j] < 0:
            t = d[i] / (d[i] - d[j])
            newV.append(vadd(V[i], vscale(t, vsub(V[j], V[i]))))
            newM.append((M[i] & M[j]) | bit)
    return newV, newM


def make_part(V, M, planes, meta=None):
    p = Part()
    p.V, p.M, p.planes, p.meta = V, M, planes, meta
    p.fV = [[fl(x) for x in v] for v in V]
    p.box = [(min(f[k] for f in p.fV), max(f[k] for f in p.fV)) for k in range(3)]
    p.edges = None
    return p


BOX_MARGIN = 1e-7  # float coordinates of the charts are below 1e3 and exact to 1e-12: a margin of 1e-7 is rigorous


def boxes_apart(a, b):
    for k in range(3):
        if a.box[k][1] < b.box[k][0] - BOX_MARGIN or b.box[k][1] < a.box[k][0] - BOX_MARGIN:
            return True
    return False


def part_edges(p):
    if p.edges is None:
        p.edges = adjacency(p.M)
    return p.edges


def _exact_separates(a, b, n):
    """Is max over a of n . y <= min over b of n . y, exactly?"""
    hi = None
    for v in a.V:
        x = vdot(n, v)
        if hi is None or x > hi:
            hi = x
    for v in b.V:
        if qs(vdot(n, v) - hi) < 0:
            return False
    return True


def parts_separated(a, b, stats=None):
    """True only with an exact separating plane between two convex parts (interiors may touch): a face normal of either, or the cross product of
    an edge of each. False when none is found (the parts overlap, or the search gave up)."""
    if boxes_apart(a, b):
        return True
    cands = []

    def float_gap(nf):
        pa = [nf[0] * f[0] + nf[1] * f[1] + nf[2] * f[2] for f in a.fV]
        pb = [nf[0] * f[0] + nf[1] * f[1] + nf[2] * f[2] for f in b.fV]
        # gap > 0: a below b along n; gap2: b below a
        return min(pb) - max(pa), min(pa) - max(pb)

    def consider(getn, nf):
        g1, g2 = float_gap(nf)
        scale = math.sqrt(nf[0] ** 2 + nf[1] ** 2 + nf[2] ** 2)
        if scale < 1e-12:
            return
        g1, g2 = g1 / scale, g2 / scale
        if g1 >= -1e-7:
            cands.append((g1, getn, True))
        if g2 >= -1e-7:
            cands.append((g2, getn, False))

    for part in (a, b):
        for k, (N, h) in enumerate(part.planes):
            if sum(1 for m in part.M if (m >> k) & 1) < 3:
                continue
            consider(lambda N=N: N, [fl(x) for x in N])
    cands.sort(key=lambda t: -t[0])
    tried = 0
    for g, getn, a_below in cands[:12]:
        tried += 1
        n = getn()
        if (_exact_separates(a, b, n) if a_below else _exact_separates(b, a, n)):
            return True
        if stats is not None:
            stats["sat_exact_fail"] = stats.get("sat_exact_fail", 0) + 1
    # edge-edge axes
    cands = []
    ea, eb = part_edges(a), part_edges(b)
    for i, j in ea:
        da = [a.fV[i][k] - a.fV[j][k] for k in range(3)]
        for p, q in eb:
            db = [b.fV[p][k] - b.fV[q][k] for k in range(3)]
            nf = [da[1] * db[2] - da[2] * db[1], da[2] * db[0] - da[0] * db[2], da[0] * db[1] - da[1] * db[0]]
            consider(lambda i=i, j=j, p=p, q=q: cross(vsub(a.V[i], a.V[j]), vsub(b.V[p], b.V[q])), nf)
    cands.sort(key=lambda t: -t[0])
    for g, getn, a_below in cands[:40]:
        n = getn()
        if (_exact_separates(a, b, n) if a_below else _exact_separates(b, a, n)):
            return True
        if stats is not None:
            stats["sat_exact_fail"] = stats.get("sat_exact_fail", 0) + 1
    return False


# ------------------------------------------------------------------ the check


class Verdict:
    def __init__(self, strict, eps):
        self.strict = strict
        self.eps = eps
        self.reasons = []
        self.stats = {"chains": 0, "parts": 0, "pairs": 0, "nodes_skipped": 0}
        self.margins = {}

    @property
    def ok(self):
        return not self.reasons

    def fail(self, msg):
        if len(self.reasons) < 50:
            self.reasons.append(msg)


def _lift_mat(m):
    return [[QM.of(x) for x in row] for row in m]


def _cell_planes(cont, ci):
    """The facet planes of cell ci in its own chart, as (N, h) with the cell on the side N . y <= h (QM)."""
    c = cont.cells[ci]
    if hasattr(c, "qplanes"):
        return c.qplanes
    out = []
    for N, q in G.facet_planes(cont, c):
        Nq = qvec(N)
        out.append((Nq, vdot(Nq, qvec(q))))
    c.qplanes = out
    return out


def _lift_vec(v):
    return [QM.of(x) for x in v]


def segment_window_interval(win, a, b):
    """Parameters [t0, t1] of a + t (b - a), t in [0, 1], on the sides N . y <= h of every plane (N, h) of win; None when empty."""
    one, zero = QM(1), QM(0)
    t0, t1 = zero, one
    for N, h in win:
        fa = vdot(N, a) - h
        fb = vdot(N, b) - h
        sa, sb = qs(fa), qs(fb)
        if sa <= 0 and sb <= 0:
            continue
        if sa > 0 and sb > 0:
            return None
        t = fa / (fa - fb)
        if sa > 0:  # entering the allowed side
            if t > t0:
                t0 = t
        else:
            if t < t1:
                t1 = t
        if t0 > t1:
            return None
    return t0, t1


def segment_meets_interior(planes, a, b, t0, t1):
    """Does the open interior {y : N . y < h for every plane} meet the segment a + t (b - a), t in [t0, t1]? Exact."""
    lo = None  # t > lo
    hi = None  # t < hi
    for N, h in planes:
        alpha = vdot(N, a) - h
        beta = vdot(N, vsub(b, a))
        sb = qs(beta)
        if sb == 0:
            if qs(alpha) >= 0:
                return False
            continue
        t = -alpha / beta
        if sb > 0:
            if hi is None or t < hi:
                hi = t
        else:
            if lo is None or t > lo:
                lo = t
        if lo is not None and hi is not None and not lo < hi:
            return False
    # (lo, hi) is the open set; meet [t0, t1]: lo < t1 and hi > t0
    if lo is not None and not lo < t1:
        return False
    if hi is not None and not hi > t0:
        return False
    return t0 <= t1


def check_poses(cont, shape, s, poses, strict=True, eps=None, max_nodes=MAX_NODES, full=False):
    """Decides validity of the exact poses (Pose objects); `eps` is the inflation factor (default EPS_STRICT when strict, else 0)."""
    if eps is None:
        eps = EPS_STRICT if strict else Fraction(0)
    v = Verdict(strict, eps)
    P = len(poses)
    e2 = cont.e2
    # squared reach in base units: ((1 + eps) circumradius e0 / s)^2
    R2 = shape.circ2 * (e2 / (s * s)) * ((1 + eps) * (1 + eps))
    parts_by_cell = {}
    nodes_total = 0
    hinge_fail = set()
    for i, pose in enumerate(poses):
        closed = getattr(pose, "closed", None) or {pose.cell: pose.lam}
        lam0_a, verts_a = piece_vertices_in_anchor(cont, shape, s, pose, eps)
        for ci, lam_c in closed.items():
            if ci == pose.cell:
                verts_run = verts_a
                lam_run = lam0_a
            else:
                fi = next((k for k, f in enumerate(cont.cells[ci].facets) if f["nbr"] == pose.cell), None)
                if fi is None:
                    v.fail(f"piece {i}: its centre lies on an edge of the polytope")
                    return v
                A, t = cont.ginv(cont.cells[ci], fi)
                Aq = _lift_mat(A)
                tq = _lift_vec(t)
                verts_run = [vadd(matvec(Aq, p), tq) for p in verts_a]
                lam_run = vadd(matvec(Aq, lam0_a), tq)
            try:
                nodes = G.enumerate_chains(cont, ci, [x if isinstance(x, Q5) else x for x in lam_c], R2, max_nodes)
            except G.ChainLimit as e:
                raise PlatonicError(f"piece {i}: {e}") from None
            nodes_total += len(nodes)
            _develop_run(cont, shape, v, i, ci, nodes, lam_run, verts_run, parts_by_cell, hinge_fail)
            if v.reasons and not full:
                return v
    v.stats["chains"] = nodes_total
    v.stats["parts"] = sum(len(x) for x in parts_by_cell.values())
    if v.reasons and not full:
        return v
    # pairs in every cell
    for cell, plist in parts_by_cell.items():
        for x in range(len(plist)):
            a = plist[x]
            for y in range(x + 1, len(plist)):
                b = plist[y]
                if a.meta[0] == b.meta[0] and a.meta[1] == b.meta[1] and a.meta[2] == b.meta[2]:
                    continue
                if boxes_apart(a, b):
                    continue
                v.stats["pairs"] += 1
                if not parts_separated(a, b, v.stats):
                    i, j = a.meta[0], b.meta[0]
                    if i == j:
                        v.fail(f"piece {i} overlaps itself in cell {cell}")
                    else:
                        v.fail(f"pieces {min(i, j)} and {max(i, j)} overlap in cell {cell}")
                    if not full:
                        return v
    return v


def _develop_run(cont, shape, v, i, ci, nodes, lam, verts, parts_by_cell, hinge_fail):
    """Parts and hinge tests of one run (one anchor cell) of piece i."""
    nf = len(shape.faces)
    for ni, node in enumerate(nodes):
        cc = cont.cells[node.cell]
        A, t = node.A, node.t
        Ai = _lift_mat(inv3(A))
        tq = _lift_vec(t)
        lamq = lam
        cen = matvec(Ai, vsub(lamq, tq))
        P = [matvec(Ai, vsub(p, tq)) for p in verts]
        # quick float prefilter on the cell's bounding box (rigorous with the margin)
        pbox = [(min(fl(p[k]) for p in P), max(fl(p[k]) for p in P)) for k in range(3)]
        cbox = getattr(cc, "fbox", None)
        if cbox is None:
            cv = [[float(x) for x in cc.chart_v[idx]] for idx in cc.ids]
            cbox = [(min(p[k] for p in cv), max(p[k] for p in cv)) for k in range(3)]
            cc.fbox = cbox
        # window planes pulled back into this cell's chart: (-A^T N) . y <= N . (t - lam)
        win = []
        for N in node.window:
            Nq = _lift_vec(N)
            AtN = [vdot([QM.of(A[r][cc_]) for r in range(3)], Nq) for cc_ in range(3)]
            win.append(([-x for x in AtN], vdot(Nq, vsub(tq, lamq))))
        # the piece's face planes (outward)
        planes, M = [], []
        for f in shape.faces:
            a, b, c = f["tri"]
            Nf = cross(vsub(P[b], P[a]), vsub(P[c], P[a]))
            h = vdot(Nf, P[a])
            if qs(vdot(Nf, cen) - h) > 0:
                Nf = [-x for x in Nf]
                h = -h
            planes.append((Nf, h))
        # (1) the hinges of this cell inside the closed window against the open interior of the piece
        in_box = all(not (pbox[k][1] < cbox[k][0] - BOX_MARGIN or cbox[k][1] < pbox[k][0] - BOX_MARGIN) for k in range(3))
        if in_box:
            for a_id, b_id in cc.hinges:
                ha = [float(x) for x in cc.chart_v[a_id]]
                hb = [float(x) for x in cc.chart_v[b_id]]
                if any(max(ha[k], hb[k]) < pbox[k][0] - BOX_MARGIN or min(ha[k], hb[k]) > pbox[k][1] + BOX_MARGIN for k in range(3)):
                    continue
                ea, eb = _lift_vec(cc.chart_v[a_id]), _lift_vec(cc.chart_v[b_id])
                iv = segment_window_interval(win, ea, eb)
                if iv is None:
                    continue
                if segment_meets_interior(planes, ea, eb, iv[0], iv[1]):
                    v.fail(f"an edge of the polytope passes through the interior of piece {i} (cell {node.cell})")
                    hinge_fail.add(i)
        # (2) the part in this cell
        if not in_box:
            v.stats["nodes_skipped"] += 1
            continue
        M = list(shape.vmask)
        pl = list(planes)
        cur = (P, M)
        for N, h in _cell_planes(cont, node.cell):
            cur = clip(cur[0], cur[1], pl, N, h)
            if cur is None:
                break
        if cur is None:
            continue
        for N, h in win:
            cur = clip(cur[0], cur[1], pl, N, h)
            if cur is None:
                break
        if cur is None:
            continue
        if not affine_rank3(cur[0]):
            continue
        part = make_part(cur[0], cur[1], pl, meta=(i, ci, ni))
        parts_by_cell.setdefault(node.cell, []).append(part)


# ------------------------------------------------------------------ the front end used by certify.py


def exact_shape_volume(name):
    """Exact volume of the unit-edge solid as (a, b, r): a + b sqrt r with r in {2, 5} (b = 0 when rational)."""
    F = Fraction
    return {
        "tetrahedron": (F(0), F(1, 12), 2),
        "cube": (F(1), F(0), 2),
        "octahedron": (F(0), F(1, 3), 2),
        "dodecahedron": (F(15, 4), F(7, 4), 5),
        "icosahedron": (F(15, 12), F(5, 12), 5),
    }[name]


def exact_bounds(polytope, piece, P):
    """The proven lower bounds on s for P pieces (BOUNDARY.md 4C.4) as (kind, text, check): check(x) decides that the rational x is not above the
    bound. Volume: P V_piece <= V_boundary s^3. One piece: s >= r_in / rho (the inscribed ball must be a valid ball)."""
    F = Fraction
    vol = {
        "tesseract": (F(8), F(0), 2),
        "5cell": (F(0), F(5, 12), 2),
        "16cell": (F(0), F(4, 3), 2),
        "24cell": (F(0), F(8), 2),
        "600cell": (F(0), F(50), 2),
        "120cell": (F(450), F(210), 5),
    }[polytope]
    rho2 = {"tesseract": Q5(F(1, 2)), "5cell": Q5(F(1, 8)), "16cell": Q5(F(1, 8)), "600cell": Q5(F(1, 8)), "24cell": Q5(F(1, 4)), "120cell": Q5(F(7, 8), F(3, 8))}[polytope]
    pa, pb, pr = exact_shape_volume(piece)
    sh = get_shape(piece)
    lo2, hi2, lo5, hi5 = G.SQRT2_LO, G.SQRT2_HI, G.SQRT5_LO, G.SQRT5_HI

    def enc(a, b, r):
        # (lo, hi) of a + b sqrt r
        l, h = (lo2, hi2) if r == 2 else (lo5, hi5)
        return a + b * l, a + b * h

    def volume_check(x):
        # x^3 V_boundary <= P V_piece, decided with the enclosures: a conservative test that says "not above" only when sure
        vb_hi = enc(*vol)[1]
        vp_lo = enc(pa, pb, pr)[0]
        return x ** 3 * vb_hi <= P * vp_lo

    def one_check(x):
        # x <= r_in / rho  <=>  x^2 rho^2 <= r_in^2
        return qs(rho2 * (x * x) - sh.in2) <= 0

    return [("volume", "(P V_piece / V_boundary)^(1/3)", volume_check), ("one piece", "r_in / rho", one_check)]


def check_pieces(polytope, piece, P, config, to_fraction, strict=True, eps=None, max_nodes=MAX_NODES, full=False):
    """Decides validity (strict: strict validity) of a float configuration of P pieces. Raises PlatonicError for malformed input or work beyond
    the caps."""
    if isinstance(P, bool) or not isinstance(P, int) or not 1 <= P <= 9999:
        raise PlatonicError("P must be a whole number from 1 to 9999")
    cfg = parse_config(polytope, piece, P, config, to_fraction)
    shape = get_shape(piece)
    cont = G.get_container(polytope)
    s = cfg["s"]
    # cheap first: the two lower bounds decide a tiny s at once. Strict: s equal to the bound is invalid.
    sx = s if strict else s * (1 + Fraction(1, 10 ** 40))
    for kind, text, below in exact_bounds(polytope, piece, P):
        if below(sx):
            v = Verdict(strict, eps if eps is not None else Fraction(0))
            v.fail(f"s is too small for {P} {piece} pieces: it must be {'more than' if strict else 'at least'} the {kind} bound {text}")
            return v
    try:
        poses = poses_from_config(cont, cfg)
    except ValueError as e:
        v = Verdict(strict, Fraction(0))
        v.fail(str(e))
        return v
    return check_poses(cont, shape, s, poses, strict=strict, eps=eps, max_nodes=max_nodes, full=full)


def summary(v):
    out = {"checker": "python/boundary/platonic.py", "version": "boundary-platonic 1", "strict": v.strict, "valid": v.ok}
    out["inflation"] = float(v.eps)
    out["chains"] = v.stats.get("chains")
    out["parts"] = v.stats.get("parts")
    out["pairs_tested"] = v.stats.get("pairs")
    return out


# ------------------------------------------------------------------ measuring the reading


def reading_discrepancy(polytope, piece, config, to_fraction):
    """The largest distance, over the piece's vertices, between the exact pose this module reads from a float configuration and the engine's own
    float computation of the same vertices (ambient coordinates at the configuration's s), divided by s. Floats on both sides; the exact side is
    rounded once to binary64."""
    cfg = parse_config(polytope, piece, len(config["pieces"]), config, to_fraction)
    cont = G.get_container(polytope)
    shape = get_shape(piece)
    poses = poses_from_config(cont, cfg)
    s = cfg["s"]
    worst = 0.0
    for pose, p in zip(poses, cfg["pieces"]):
        c = cont.cells[pose.cell]
        lam0, verts = piece_vertices_in_anchor(cont, shape, s, pose, Fraction(0))
        fv = cont.fvertices
        Re = quat_matrix_float([float(x) for x in p["q"]])
        for vs, w in zip(shape.V, verts):
            lam = [fl(x) for x in w]
            ex = [0.0] * 4
            for i in range(4):
                ex[i] = fv[c.O_idx][i] + sum(lam[j] * (fv[c.tri[j]][i] - fv[c.O_idx][i]) for j in range(3))
            # engine: s * origin + sum_k frame_k (u_k + (R_e v)_k)
            vstd = [fl(x) for x in vs]
            loc = [float(p["u"][k]) + sum(Re[k][j] * vstd[j] for j in range(3)) for k in range(3)]
            en = [float(s) * c.forigin[i] + sum(loc[k] * c.fframe[k][i] for k in range(3)) for i in range(4)]
            worst = max(worst, math.sqrt(sum((float(s) * ex[i] - en[i]) ** 2 for i in range(4))) / float(s))
    return worst


def main(argv):
    import json
    import sys

    if len(argv) != 2:
        print("usage: python3 python/boundary/platonic.py boundary/<polytope>/<piece>/<P> <config.json>")
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
