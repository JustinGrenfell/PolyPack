/**
 * The Platonic solids as pieces of the boundary space (3D on 4D), unit edge, centred at the origin.
 *
 * A piece's standard pose is the solid below; a configuration's `q` rotates it about its centre in
 * the anchor cell's frame. The cube is the one that was there first (parts.ts CUBE_VERTICES); the
 * others follow the same contract: vertex coordinates exact in Q(sqrt 2) (tetrahedron, octahedron)
 * or Q(sqrt 5) (dodecahedron, icosahedron), faces counter-clockwise seen from outside.
 *
 *  - tetrahedron: (1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1) / (2 sqrt 2);
 *  - cube: (+-1/2)^3;
 *  - octahedron: +-e_k / sqrt 2;
 *  - dodecahedron: (+-phi/2)^3 and the cyclic shifts of (0, +-1/2, +-phi^2/2);
 *  - icosahedron: the cyclic shifts of (0, +-1/2, +-phi/2);
 * phi = (1 + sqrt 5) / 2. python/boundary/platonic.py repeats these coordinates in exact arithmetic.
 */
import { V3, cross3, dot3, sub3, unit3 } from "./linalg";
import type { PolytopePieceName } from "./types";

const PHI = (1 + Math.sqrt(5)) / 2;

export interface PieceShape {
  name: PolytopePieceName;
  vertices: V3[];
  /** Faces as vertex indices, counter-clockwise seen from outside. */
  faces: number[][];
  /** Outward unit normals of the faces. */
  normals: V3[];
  /** Distance from the centre to every face. */
  inradius: number;
  /** Distance from the centre to every vertex. */
  circumradius: number;
  /** Distance from the centre to every edge. */
  midradius: number;
  volume: number;
  /** An orientation-reversing symmetry of the solid (a matrix, row-major): R S = a rotation image of R. */
  mirror: number[];
  /** Edges as vertex index pairs. */
  edges: [number, number][];
  /** One unit direction per class of parallel edges (the separating-axis test's edge axes). */
  edgeDirs: V3[];
}

/** Faces of the convex hull of `verts` (counter-clockwise from outside), found by planes through triples. */
function hullFaces(verts: readonly V3[]): { faces: number[][]; normals: V3[] } {
  const n = verts.length;
  const found: { normal: V3; ids: number[] }[] = [];
  for (let a = 0; a < n; a++)
    for (let b = a + 1; b < n; b++)
      for (let c = b + 1; c < n; c++) {
        const nr = cross3(sub3(verts[b], verts[a]), sub3(verts[c], verts[a]));
        const len = Math.hypot(...nr);
        if (len < 1e-9) continue;
        const nu = unit3(nr);
        const h = dot3(nu, verts[a]);
        let pos = 0;
        let neg = 0;
        for (let i = 0; i < n; i++) {
          const d = dot3(nu, verts[i]) - h;
          if (d > 1e-9) pos++;
          else if (d < -1e-9) neg++;
        }
        if (pos && neg) continue;
        const out: V3 = pos ? [-nu[0], -nu[1], -nu[2]] : nu;
        if (found.some((f) => dot3(f.normal, out) > 1 - 1e-9)) continue;
        const hh = dot3(out, verts[a]);
        const ids = verts.map((v, i) => i).filter((i) => Math.abs(dot3(out, verts[i]) - hh) < 1e-9);
        found.push({ normal: out, ids });
      }
  const faces = found.map((f) => {
    // order counter-clockwise around the normal
    const c: V3 = [0, 0, 0];
    for (const i of f.ids) for (let k = 0; k < 3; k++) c[k] += verts[i][k] / f.ids.length;
    const e1 = unit3(sub3(verts[f.ids[0]], c));
    const e2 = cross3(f.normal, e1);
    const ang = (i: number) => {
      const d = sub3(verts[i], c);
      return Math.atan2(dot3(d, e2), dot3(d, e1));
    };
    return f.ids.slice().sort((x, y) => ang(x) - ang(y));
  });
  return { faces, normals: found.map((f) => f.normal) };
}

function edgesOf(faces: number[][]): [number, number][] {
  const set = new Set<string>();
  const out: [number, number][] = [];
  for (const f of faces)
    for (let i = 0; i < f.length; i++) {
      const a = Math.min(f[i], f[(i + 1) % f.length]);
      const b = Math.max(f[i], f[(i + 1) % f.length]);
      const k = `${a},${b}`;
      if (!set.has(k)) {
        set.add(k);
        out.push([a, b]);
      }
    }
  return out.sort((p, q) => p[0] - q[0] || p[1] - q[1]);
}

function signs(base: V3): V3[] {
  const out: V3[] = [];
  for (let m = 0; m < 8; m++) {
    const v: V3 = [m & 1 ? -base[0] : base[0], m & 2 ? -base[1] : base[1], m & 4 ? -base[2] : base[2]];
    if (!out.some((w) => w[0] === v[0] && w[1] === v[1] && w[2] === v[2])) out.push(v);
  }
  return out;
}

function cyclic(v: V3): V3[] {
  return [v, [v[1], v[2], v[0]], [v[2], v[0], v[1]]];
}

function buildVertices(name: PolytopePieceName): V3[] {
  switch (name) {
    case "tetrahedron": {
      const k = 1 / (2 * Math.SQRT2);
      return ([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]] as V3[]).map((v) => [v[0] * k, v[1] * k, v[2] * k] as V3);
    }
    case "cube":
      return Array.from({ length: 8 }, (_, v) => [v & 1 ? 0.5 : -0.5, v & 2 ? 0.5 : -0.5, v & 4 ? 0.5 : -0.5] as V3);
    case "octahedron": {
      const r = Math.SQRT1_2;
      return [[r, 0, 0], [-r, 0, 0], [0, r, 0], [0, -r, 0], [0, 0, r], [0, 0, -r]];
    }
    case "dodecahedron": {
      const out: V3[] = signs([PHI / 2, PHI / 2, PHI / 2]);
      for (const c of cyclic([0, 0.5, (PHI * PHI) / 2])) for (const v of signs(c)) out.push(v);
      return out;
    }
    case "icosahedron": {
      const out: V3[] = [];
      for (const c of cyclic([0, 0.5, PHI / 2])) for (const v of signs(c)) out.push(v);
      return out;
    }
  }
}

const MIRROR_SWAP_XY = [0, 1, 0, 1, 0, 0, 0, 0, 1];
const MIRROR_INVERSION = [-1, 0, 0, 0, -1, 0, 0, 0, -1];

function build(name: PolytopePieceName): PieceShape {
  const vertices = buildVertices(name);
  // the cube keeps its published face order (parts.ts CUBE_FACES: -x, +x, -y, +y, -z, +z)
  const hull =
    name === "cube"
      ? {
          faces: [[0, 4, 6, 2], [1, 3, 7, 5], [0, 1, 5, 4], [2, 6, 7, 3], [0, 2, 3, 1], [4, 5, 7, 6]],
          normals: [[-1, 0, 0], [1, 0, 0], [0, -1, 0], [0, 1, 0], [0, 0, -1], [0, 0, 1]] as V3[],
        }
      : hullFaces(vertices);
  const edges = edgesOf(hull.faces);
  const edgeDirs: V3[] = [];
  for (const [a, b] of edges) {
    const d = unit3(sub3(vertices[b], vertices[a]));
    if (!edgeDirs.some((e) => Math.hypot(...cross3(e, d)) < 1e-9)) edgeDirs.push(d);
  }
  const inradius = dot3(hull.normals[0], vertices[hull.faces[0][0]]);
  const circumradius = name === "cube" ? Math.sqrt(3) / 2 : Math.hypot(...vertices[0]);
  const [ea, eb] = edges[0];
  const mid: V3 = [(vertices[ea][0] + vertices[eb][0]) / 2, (vertices[ea][1] + vertices[eb][1]) / 2, (vertices[ea][2] + vertices[eb][2]) / 2];
  // volume: sum over faces of (1/3) area * inradius
  let vol = 0;
  hull.faces.forEach((f, k) => {
    let a = 0;
    for (let i = 1; i + 1 < f.length; i++) a += Math.hypot(...cross3(sub3(vertices[f[i]], vertices[f[0]]), sub3(vertices[f[i + 1]], vertices[f[0]]))) / 2;
    vol += (a * dot3(hull.normals[k], vertices[f[0]])) / 3;
  });
  return {
    name,
    vertices,
    faces: hull.faces,
    normals: hull.normals,
    inradius,
    circumradius,
    midradius: Math.hypot(...mid),
    volume: vol,
    mirror: name === "tetrahedron" ? MIRROR_SWAP_XY : MIRROR_INVERSION,
    edges,
    edgeDirs,
  };
}

export const POLYTOPE_PIECES: readonly PolytopePieceName[] = ["tetrahedron", "cube", "octahedron", "dodecahedron", "icosahedron"];

const SHAPES = new Map<PolytopePieceName, PieceShape>();

/** The shape of a Platonic piece (built once). */
export function pieceShape(name: PolytopePieceName): PieceShape {
  let s = SHAPES.get(name);
  if (!s) {
    s = build(name);
    SHAPES.set(name, s);
  }
  return s;
}

export function isPolytopePiece(x: unknown): x is PolytopePieceName {
  return typeof x === "string" && (POLYTOPE_PIECES as readonly string[]).includes(x);
}

function inv3m(m: readonly number[]): number[] | null {
  const [a, b, c, d, e, f, g, h, i] = m;
  const det = a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g);
  if (Math.abs(det) < 1e-12) return null;
  return [
    (e * i - f * h) / det, (c * h - b * i) / det, (b * f - c * e) / det,
    (f * g - d * i) / det, (a * i - c * g) / det, (c * d - a * f) / det,
    (d * h - e * g) / det, (b * g - a * h) / det, (a * e - b * d) / det,
  ];
}

/**
 * A proper rotation matrix (row-major) taking the standard `shape` onto the solid with the given vertices (the same solid at unit
 * edge, in any pose, centre anywhere), or null when there is none: used to place a piece exactly onto a cell of its own shape.
 */
export function rotationOnto(shape: PieceShape, vertices: readonly V3[]): number[] | null {
  const n = vertices.length;
  if (n !== shape.vertices.length) return null;
  const c: V3 = [0, 0, 0];
  for (const v of vertices) for (let k = 0; k < 3; k++) c[k] += v[k] / n;
  const T = vertices.map((v) => sub3(v, c));
  const S = shape.vertices;
  // three vertices of the standard solid that are linearly independent
  let tri: [number, number, number] | null = null;
  let Sinv: number[] | null = null;
  outer: for (let a = 0; a < S.length; a++)
    for (let b = a + 1; b < S.length; b++)
      for (let d = b + 1; d < S.length; d++) {
        const inv = inv3m([S[a][0], S[b][0], S[d][0], S[a][1], S[b][1], S[d][1], S[a][2], S[b][2], S[d][2]]);
        if (inv) {
          tri = [a, b, d];
          Sinv = inv;
          break outer;
        }
      }
  if (!tri || !Sinv) return null;
  const dist = (p: V3, q: V3) => Math.hypot(p[0] - q[0], p[1] - q[1], p[2] - q[2]);
  for (let j0 = 0; j0 < n; j0++)
    for (let j1 = 0; j1 < n; j1++) {
      if (j1 === j0 || Math.abs(dist(T[j0], T[j1]) - dist(S[tri[0]], S[tri[1]])) > 1e-9) continue;
      for (let j2 = 0; j2 < n; j2++) {
        if (j2 === j0 || j2 === j1) continue;
        if (Math.abs(dist(T[j0], T[j2]) - dist(S[tri[0]], S[tri[2]])) > 1e-9 || Math.abs(dist(T[j1], T[j2]) - dist(S[tri[1]], S[tri[2]])) > 1e-9) continue;
        // R S3 = T3  =>  R = T3 S3^-1 (columns)
        const T3 = [T[j0][0], T[j1][0], T[j2][0], T[j0][1], T[j1][1], T[j2][1], T[j0][2], T[j1][2], T[j2][2]];
        const R: number[] = new Array(9).fill(0);
        for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) for (let k = 0; k < 3; k++) R[3 * i + j] += T3[3 * i + k] * Sinv[3 * k + j];
        const det = R[0] * (R[4] * R[8] - R[5] * R[7]) - R[1] * (R[3] * R[8] - R[5] * R[6]) + R[2] * (R[3] * R[7] - R[4] * R[6]);
        if (det < 0) continue;
        // every vertex of the standard solid must land on a vertex of the target
        const ok = S.every((v) => {
          const w: V3 = [R[0] * v[0] + R[1] * v[1] + R[2] * v[2], R[3] * v[0] + R[4] * v[1] + R[5] * v[2], R[6] * v[0] + R[7] * v[1] + R[8] * v[2]];
          return T.some((t) => dist(t, w) < 1e-9);
        });
        if (ok) return R;
      }
    }
  return null;
}
