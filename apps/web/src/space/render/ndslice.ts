/**
 * The geometry behind the board's drawings of the nD in nD space (docs/spaces/NDIM.md): no DOM, no three.js.
 *
 *  - the placed pieces and container of a configuration (world vertices, or a centre and a radius for a ball);
 *  - the slice by the hyperplane x_k = c: a convex polytope cut by a hyperplane is a convex polytope of one dimension
 *    less, whose vertices are the polytope's vertices on the hyperplane and the points where its edges cross it; a ball
 *    of radius r about a centre with x_k = a cuts to a ball of radius sqrt(r^2 - (c - a)^2) about the centre with the
 *    k-th coordinate dropped;
 *  - the orthographic projection onto two axes (the 2D hull of the projected vertices);
 *  - the short texts of a table row: a centre and how far a piece is turned.
 *
 * A polytope's edges are not stored in the catalogue (they are the vertex pairs at distance 1), so they are found once
 * per polytope here.
 */
import { BALL, type NdConfig, type NdProblem, type Polytope, getPolytope, loadPolytope, readPose, placeVertices } from "@polypack/ndim";

/** A ball piece or container: a centre and a radius. */
export interface PlacedBall {
  kind: "ball";
  centre: Float64Array;
  radius: number;
}

/** A polytope piece or container placed in the world: its vertices (nV x d, flat), and the rotation of a piece. */
export interface PlacedPolytope {
  kind: "polytope";
  d: number;
  V: Float64Array;
  /** Index pairs of the unit edges (scaled with the polytope). */
  edges: readonly (readonly [number, number])[];
  centre: Float64Array;
  /** Row-major d x d rotation; the identity for a container. */
  R: Float64Array;
}

export type Placed = PlacedBall | PlacedPolytope;

/** What a problem's shapes need before a configuration can be placed: the loaded polytopes. */
export async function loadShapes(p: Pick<NdProblem, "piece" | "container">): Promise<void> {
  for (const c of new Set([p.piece, p.container])) if (c !== BALL) await loadPolytope(c);
}

const edgeCache = new Map<string, readonly (readonly [number, number])[]>();

/** The unit edges of a loaded polytope: the vertex pairs at distance 1 (to 1e-6; the catalogue's edges are all exactly 1). */
export function edgesOf(poly: Polytope): readonly (readonly [number, number])[] {
  const hit = edgeCache.get(poly.code);
  if (hit) return hit;
  const { V, nV, d } = poly;
  const out: [number, number][] = [];
  for (let a = 0; a < nV; a++)
    for (let b = a + 1; b < nV; b++) {
      let s = 0;
      for (let i = 0; i < d; i++) {
        const x = V[a * d + i] - V[b * d + i];
        s += x * x;
      }
      if (Math.abs(s - 1) < 1e-6) out.push([a, b]);
    }
  edgeCache.set(poly.code, out);
  return out;
}

/** The container of a configuration: s times the polytope about the origin, or the ball of diameter s. */
export function placeContainer(p: NdProblem, s: number): Placed {
  const centre = new Float64Array(p.d);
  if (p.container === BALL) return { kind: "ball", centre, radius: s / 2 };
  const poly = getPolytope(p.container);
  const R = new Float64Array(p.d * p.d);
  for (let i = 0; i < p.d; i++) R[i * p.d + i] = 1;
  return { kind: "polytope", d: p.d, V: poly.V.map((x) => x * s), edges: edgesOf(poly), centre, R };
}

/** The pieces of a configuration (the polytopes must be loaded). A piece that does not read is skipped and reported. */
export function placePieces(p: NdProblem, config: NdConfig): { pieces: (Placed | null)[]; errors: string[] } {
  const errors: string[] = [];
  const isBall = p.piece === BALL;
  const poly = isBall ? null : getPolytope(p.piece);
  const edges = poly ? edgesOf(poly) : [];
  const pieces = config.pieces.map((pose, i): Placed | null => {
    const r = readPose(p.d, pose, isBall);
    if (typeof r === "string") {
      errors.push(`piece ${i + 1}: ${r}`);
      return null;
    }
    if (isBall || !poly) return { kind: "ball", centre: r.t, radius: 0.5 };
    return {
      kind: "polytope",
      d: p.d,
      V: placeVertices(poly.V, poly.nV, p.d, r.t, r.R!),
      edges,
      centre: r.t,
      R: r.R!,
    };
  });
  return { pieces, errors };
}

/** The range of x_k over a container: [lo, hi]. */
export function extentAlong(c: Placed, k: number): [number, number] {
  if (c.kind === "ball") return [c.centre[k] - c.radius, c.centre[k] + c.radius];
  const d = c.d;
  let lo = Infinity;
  let hi = -Infinity;
  for (let v = 0; v < c.V.length / d; v++) {
    const x = c.V[v * d + k];
    if (x < lo) lo = x;
    if (x > hi) hi = x;
  }
  return [lo, hi];
}

/** A slice of a ball: a ball one dimension down, or null when the hyperplane misses it (a touch is empty too). */
export interface BallSlice {
  kind: "ball";
  centre: number[];
  radius: number;
}

/** A slice of a polytope: the points whose hull is the slice (one dimension down). */
export interface PolytopeSlice {
  kind: "polytope";
  points: number[][];
}

export type Slice = BallSlice | PolytopeSlice;

const without = (x: ArrayLike<number>, k: number): number[] => Array.from(x).filter((_, i) => i !== k);

/** Distance below which a vertex counts as lying on the hyperplane, and two cut points as one. */
export const SLICE_TOL = 1e-9;

/**
 * The slice of a placed shape by x_k = c, or null when it is empty. A polytope that only touches the hyperplane (at a
 * vertex or along a face) has a slice of its touching points, which are listed; a ball that only touches is empty.
 */
export function sliceAt(shape: Placed, k: number, c: number): Slice | null {
  if (shape.kind === "ball") {
    const a = shape.centre[k];
    const r2 = shape.radius * shape.radius - (c - a) * (c - a);
    if (!(r2 > SLICE_TOL * SLICE_TOL)) return null;
    return { kind: "ball", centre: without(shape.centre, k), radius: Math.sqrt(r2) };
  }
  const { d, V, edges } = shape;
  const n = V.length / d;
  const points: number[][] = [];
  const put = (x: number[]) => {
    for (const q of points) {
      let s = 0;
      for (let i = 0; i < x.length; i++) s = Math.max(s, Math.abs(q[i] - x[i]));
      if (s < SLICE_TOL) return;
    }
    points.push(x);
  };
  const row = (v: number) => V.subarray(v * d, (v + 1) * d);
  for (let v = 0; v < n; v++) if (Math.abs(V[v * d + k] - c) <= SLICE_TOL) put(without(row(v), k));
  for (const [a, b] of edges) {
    const ha = V[a * d + k] - c;
    const hb = V[b * d + k] - c;
    // an end on the hyperplane was taken as a vertex above
    if (Math.abs(ha) <= SLICE_TOL || Math.abs(hb) <= SLICE_TOL || ha * hb > 0) continue;
    const t = ha / (ha - hb);
    const x = new Array<number>(d);
    for (let i = 0; i < d; i++) x[i] = V[a * d + i] + t * (V[b * d + i] - V[a * d + i]);
    put(without(x, k));
  }
  return points.length ? { kind: "polytope", points } : null;
}

// ---------- the projection onto two axes ----------

/** The 2D hull of points, counter-clockwise (Andrew's monotone chain); collinear points are dropped. */
export function hull2d(points: readonly (readonly [number, number])[]): [number, number][] {
  const p = points.map((q) => [q[0], q[1]] as [number, number]).sort((a, b) => a[0] - b[0] || a[1] - b[1]);
  if (p.length < 3) return p;
  const cross = (o: number[], a: number[], b: number[]) => (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0]);
  const build = (list: [number, number][]) => {
    const h: [number, number][] = [];
    for (const q of list) {
      while (h.length >= 2 && cross(h[h.length - 2], h[h.length - 1], q) <= 1e-12) h.pop();
      h.push(q);
    }
    h.pop();
    return h;
  };
  const lower = build(p);
  const upper = build([...p].reverse());
  return [...lower, ...upper];
}

/** A shape projected orthographically onto axes i and j (0-based): a disk, or the hull of the projected vertices. */
export type Projection = { kind: "ball"; cx: number; cy: number; r: number } | { kind: "polygon"; points: [number, number][] };

export function projectOnto(shape: Placed, i: number, j: number): Projection {
  if (shape.kind === "ball") return { kind: "ball", cx: shape.centre[i], cy: shape.centre[j], r: shape.radius };
  const { d, V } = shape;
  const pts: [number, number][] = [];
  for (let v = 0; v < V.length / d; v++) pts.push([V[v * d + i], V[v * d + j]]);
  return { kind: "polygon", points: hull2d(pts) };
}

// ---------- the table's texts ----------

/** "0.5000" with a true minus sign, and no "-0.0000". */
export function coordText(x: number, digits = 4): string {
  const t = x.toFixed(digits);
  const z = Number(t) === 0 ? (0).toFixed(digits) : t;
  return z.replace(/^-/, "−");
}

/** A piece's centre as "(0.5000, −0.2500, …)". */
export const centreText = (c: ArrayLike<number>, digits = 4): string =>
  `(${Array.from(c, (x) => coordText(x, digits)).join(", ")})`;

/**
 * How far a polytope piece is turned from the catalogue's orientation, in a few words: "not turned" for the identity,
 * "axes permuted" for a signed permutation matrix (the pieces of a grid), else "turned about 37°". The angle is the
 * mean turn of a vector, arccos(tr R / d): the average of cos(angle(x, R x)) over the unit sphere is tr R / d.
 */
export function rotationText(R: ArrayLike<number>, d: number): string {
  let maxDev = 0;
  let permutation = true;
  let tr = 0;
  for (let i = 0; i < d; i++) {
    tr += R[i * d + i];
    let big = 0;
    for (let j = 0; j < d; j++) {
      const x = R[i * d + j];
      maxDev = Math.max(maxDev, Math.abs(x - (i === j ? 1 : 0)));
      const a = Math.abs(x);
      if (a > 1e-6 && Math.abs(a - 1) > 1e-6) permutation = false;
      if (a > 0.5) big++;
    }
    if (big !== 1) permutation = false;
  }
  if (maxDev < 1e-9) return "not turned";
  if (permutation) return "axes permuted";
  const deg = (Math.acos(Math.max(-1, Math.min(1, tr / d))) * 180) / Math.PI;
  return `turned about ${Math.round(deg)}°`;
}
