/**
 * Half-chains, windows, joins and hinge distances on the boundary of the simplex, orthoplex and demicube (d = 5..10); the float
 * engine of docs/spaces/BOUNDARY.md section 4C.2 to 4C.4. The exact twin is python/boundary/ndpoly.py.
 *
 * Development. The cell c_0 (the anchor) holds the source point u. A chain is a sequence of cells c_0, c_1, ..., c_k, each glued to the
 * previous one across a ridge. The cell c_(j+1) is unfolded into the hyperplane of c_j by the reflection
 *     sigma(x) = x - 2 (m . (x - q)) / (m . m) m,    m = n_(c_j) + n_(c_(j+1)),   q a point of the ridge,
 * the reflection about the hyperplane through the ridge with normal n1 + n2 (n1, n2 the cells' outward unit normals): it is an
 * involution of the ambient space, fixes the ridge pointwise, maps H_(c_(j+1)) onto H_(c_j) and, because the two cells are on opposite
 * sides of the ridge's (d-2)-flat, lays c_(j+1) on the far side of the ridge (the other reflection through the ridge, with normal
 * n1 - n2, folds it back over c_j). The development of c_k is T_k = sigma_1 o sigma_2 o ... (composed on the right) : actual cell c_k ->
 * the anchor's hyperplane; a node stores its linear part L and translation t (T x = L x + t) and a = T^(-1) u, the source as seen in the
 * actual hyperplane of its own cell; a child's a is sigma(a) (an involution).
 *
 * Pruning (sound: a branch is cut only when no straight path of length <= reach can follow it). The next ridge, of neighbour g, is
 * crossed by a ray from u only if a is strictly inside the half-space of g (S_g = s h_g - n_g . a > 0), and then the path is at least
 * as long as the distance from a to the ridge, which is at least
 *     sqrt( S_g^2 / |nu_g|^2 + max_k (viol_k / |tn_k|)^2 ),
 * the distance within H_c from a to the flat of the ridge plus the distance, within that flat, from the foot of the perpendicular to the
 * half-space of the cell's other facet k (viol_k > 0 how far the foot is outside it; |tn_k| the length of that facet's normal tangential
 * to the flat). Any single facet gives a lower bound on the distance to the ridge, so the maximum does too.
 *
 * Windows. A ray from u with direction e (in the anchor's hyperplane) crosses the ridges of a chain in order, each inside the closed ridge,
 * exactly when e is in the window of the chain, a polyhedral cone: for each step (cell p, crossed neighbour f, S_k = s h_k - n_k . a_p),
 * with z the direction in p's frame, n_f . z >= 0 and S_f (n_k . z) - S_k (n_f . z) <= 0 for every other neighbour k of p
 * (derived in 4C.2: the ray hits the plane of f at tau = S_f / (n_f . z) and satisfies the facet k there). Directions are carried from
 * the end cell's frame back along the chain by the reflections (z_(j) = R_j z_(j+1)), so no window is stored.
 *
 * Joins. Two half-chains, one from each source, ending in the same cell c, with sources a1, a2 as seen in c's hyperplane, give the
 * straight path of length |a1 - a2| from the first source to the second exactly when the direction w = a2 - a1 is in the window of the
 * first and beyond its last ridge (n_f . z >= S_f), and -w is in the window of the second and beyond its last ridge. A path of length
 * L has a midpoint at distance L/2 from each source, so it is found from half-chains of reach L/2.
 *
 * The step budget: every node created, every join candidate and every projection counts one step; a check that exceeds its budget
 * throws NDLimitError and is refused, never accepted.
 */
import type { NDContainer, NDCell } from "./container";

export class NDLimitError extends Error {}

export interface NDBudget {
  steps: number;
  max: number;
}
export const newNDBudget = (max: number): NDBudget => ({ steps: 0, max });
function spend(b: NDBudget, k = 1): void {
  b.steps += k;
  if (b.steps > b.max) throw new NDLimitError(`the check needs more than ${b.max} search steps`);
}

export interface NDNode {
  cell: number;
  /** Index of the parent node (-1 for the anchor) and the index in the parent's neighbour list of the ridge crossed. */
  parent: number;
  via: number;
  depth: number;
  /** The source as seen in this cell's hyperplane. */
  a: Float64Array;
  /** The development T x = L x + t of this cell into the anchor's hyperplane (L row-major D x D). */
  L: Float64Array;
  t: Float64Array;
}

const dot = (a: ArrayLike<number>, ao: number, b: ArrayLike<number>, bo: number, D: number): number => {
  let s = 0;
  for (let j = 0; j < D; j++) s += a[ao + j] * b[bo + j];
  return s;
};

/** Reflects v about the direction m (a vector of squared length mm): v - 2 (m . v) / mm m, in place. */
function reflectDir(v: Float64Array, m: Float64Array, mm: number, D: number): void {
  const c = (2 * dot(m, 0, v, 0, D)) / mm;
  for (let j = 0; j < D; j++) v[j] -= c * m[j];
}

/** m = n_c + n_g for the ridge between cell cl and its neighbour k; returns mm = m . m. */
function mirror(cl: NDCell, k: number, D: number, m: Float64Array): number {
  let mm = 0;
  for (let j = 0; j < D; j++) {
    m[j] = cl.n[j] + cl.nn[k * D + j];
    mm += m[j] * m[j];
  }
  return mm;
}

/** Relative slack below which a ridge counts as not strictly ahead of the source. */
const AHEAD = 1e-15;

/**
 * The half-chains from the point x0 (ambient, in the hyperplane of cell c0, edge s) of reach `reach`: every chain a straight path of
 * length <= reach can follow. The anchor is node 0.
 */
export function develop(
  C: NDContainer,
  s: number,
  c0: number,
  x0: Float64Array,
  reach: number,
  budget: NDBudget,
): NDNode[] {
  const D = C.D;
  const reach2 = reach * reach * (1 + 1e-12);
  const I = new Float64Array(D * D);
  for (let j = 0; j < D; j++) I[j * D + j] = 1;
  const nodes: NDNode[] = [
    { cell: c0, parent: -1, via: -1, depth: 0, a: Float64Array.from(x0), L: I, t: new Float64Array(D) },
  ];
  const m = new Float64Array(D);
  const nu = new Float64Array(D);
  const Lm = new Float64Array(D);
  spend(budget);
  for (let idx = 0; idx < nodes.length; idx++) {
    const nd = nodes[idx];
    const cl = C.cell(nd.cell);
    const F = cl.nbr.length;
    const a = nd.a;
    const parentCell = nd.parent >= 0 ? nodes[nd.parent].cell : -1;
    // S_k = s h_k - n_k . a for every neighbour, computed once per node
    const S = new Float64Array(F);
    for (let k = 0; k < F; k++) S[k] = s * cl.hn[k] - dot(cl.nn, k * D, a, 0, D);
    for (let k = 0; k < F; k++) {
      if (cl.nbr[k] === parentCell) continue;
      const Sk = S[k];
      if (!(Sk > AHEAD * Math.max(1, s))) continue;
      const nu2 = cl.nu2[k];
      if (!(nu2 > 1e-12)) continue;
      const d2 = (Sk * Sk) / nu2;
      if (d2 > reach2) continue;
      // nu = tangential part of n_k; q_l = n_l . nu for the other facets
      const nk = cl.nn;
      const c1 = dot(nk, k * D, cl.n, 0, D);
      for (let j = 0; j < D; j++) nu[j] = nk[k * D + j] - c1 * cl.n[j];
      const f = Sk / nu2;
      let extra2 = 0;
      for (let l = 0; l < F; l++) {
        if (l === k) continue;
        const q = dot(nk, l * D, nu, 0, D);
        const viol = -S[l] + f * q;
        if (viol > 0) {
          const tn2 = cl.nu2[l] - (q * q) / nu2;
          if (tn2 > 1e-12) {
            const e2 = (viol * viol) / tn2;
            if (e2 > extra2) extra2 = e2;
          }
        }
      }
      if (d2 + extra2 > reach2) continue;
      spend(budget);
      // the child: unfold neighbour k's cell into H_c
      const mm = mirror(cl, k, D, m);
      const qv = new Float64Array(D);
      for (let j = 0; j < D; j++) qv[j] = s * cl.rp[k * D + j];
      const mq = dot(m, 0, qv, 0, D);
      const a2 = new Float64Array(D);
      const ma = (2 * (dot(m, 0, a, 0, D) - mq)) / mm;
      for (let j = 0; j < D; j++) a2[j] = a[j] - ma * m[j];
      // L' = L - 2 (L m) m^T / mm, t' = t + 2 (m . q) / mm (L m)
      for (let r = 0; r < D; r++) Lm[r] = dot(nd.L, r * D, m, 0, D);
      const L2 = new Float64Array(D * D);
      for (let r = 0; r < D; r++) {
        const c = (2 * Lm[r]) / mm;
        for (let j = 0; j < D; j++) L2[r * D + j] = nd.L[r * D + j] - c * m[j];
      }
      const t2 = new Float64Array(D);
      const ct = (2 * mq) / mm;
      for (let r = 0; r < D; r++) t2[r] = nd.t[r] + ct * Lm[r];
      nodes.push({ cell: cl.nbr[k], parent: idx, via: k, depth: nd.depth + 1, a: a2, L: L2, t: t2 });
    }
  }
  return nodes;
}

/**
 * Is the direction w (in the frame of the end cell of node idx) in the closed window of the chain ending at idx, and, when `beyond`,
 * beyond its last ridge? tol is relative (to |w| times the slacks involved).
 */
export function inWindow(
  C: NDContainer,
  s: number,
  nodes: readonly NDNode[],
  idx: number,
  w: Float64Array,
  beyond: boolean,
  tol: number,
): boolean {
  const D = C.D;
  const z = Float64Array.from(w);
  const m = new Float64Array(D);
  let nw = 0;
  for (let j = 0; j < D; j++) nw += w[j] * w[j];
  nw = Math.sqrt(nw);
  let cur = idx;
  let first = true;
  while (nodes[cur].parent >= 0) {
    const p = nodes[nodes[cur].parent];
    const pc = C.cell(p.cell);
    const f = nodes[cur].via;
    // z in p's frame: reflect by the ridge's mirror
    const mm = mirror(pc, f, D, m);
    reflectDir(z, m, mm, D);
    const F = pc.nbr.length;
    const Sf = s * pc.hn[f] - dot(pc.nn, f * D, p.a, 0, D);
    const nfz = dot(pc.nn, f * D, z, 0, D);
    if (nfz < -tol * nw) return false;
    if (first && beyond && nfz < Sf - tol * nw) return false;
    first = false;
    for (let k = 0; k < F; k++) {
      if (k === f) continue;
      const Sk = s * pc.hn[k] - dot(pc.nn, k * D, p.a, 0, D);
      const lhs = Sf * dot(pc.nn, k * D, z, 0, D) - Sk * nfz;
      if (lhs > tol * nw * (Math.abs(Sf) + Math.abs(Sk))) return false;
    }
    cur = nodes[cur].parent;
  }
  return true;
}

export interface NDJoin {
  dist: number;
  a: number;
  b: number;
}

/**
 * The shortest straight path from the source of `A` to the source of `B` found by joining half-chains, if its length is below `limit`
 * (Infinity otherwise). With `same` the two sets are one ball's: the pairs of distinct half-chains (a loop). `C` and `s` as in develop.
 */
export function joinDistance(
  C: NDContainer,
  s: number,
  A: readonly NDNode[],
  B: readonly NDNode[],
  limit: number,
  budget: NDBudget,
  same: boolean,
  tol = 1e-12,
): NDJoin {
  const D = C.D;
  const byCell = new Map<number, number[]>();
  B.forEach((nd, i) => {
    let l = byCell.get(nd.cell);
    if (!l) byCell.set(nd.cell, (l = []));
    l.push(i);
  });
  const w = new Float64Array(D);
  const wn = new Float64Array(D);
  let best: NDJoin = { dist: Infinity, a: -1, b: -1 };
  const lim2 = limit * limit;
  const zero2 = (1e-10 * Math.max(1, s)) ** 2;
  for (let i = 0; i < A.length; i++) {
    const l = byCell.get(A[i].cell);
    if (!l) continue;
    for (const j of l) {
      if (same && j <= i) continue;
      let len2 = 0;
      for (let r = 0; r < D; r++) {
        w[r] = B[j].a[r] - A[i].a[r];
        len2 += w[r] * w[r];
      }
      if (!(len2 < lim2) || len2 < zero2 || len2 >= best.dist * best.dist) continue;
      spend(budget);
      if (!inWindow(C, s, A, i, w, true, tol)) continue;
      for (let r = 0; r < D; r++) wn[r] = -w[r];
      if (!inWindow(C, s, B, j, wn, true, tol)) continue;
      best = { dist: Math.sqrt(len2), a: i, b: j };
    }
  }
  return best;
}

// ---------------------------------------------------------------------------------------------------------- hinge distance

interface Row {
  v: Float64Array;
  b: number;
}

/**
 * The squared distance from a to the polytope { x : v_i . x = b_i (i in eqs), v_j . x <= b_j (j in rows) }, or Infinity when it is empty.
 * Rows are unit vectors. The projection onto the flat of the equalities decides: when it satisfies every inequality it is the nearest
 * point; otherwise the nearest point lies on the boundary of at least one violated inequality (it is the nearest point of the
 * polytope to the projection, whose normal-cone vector has a positive inner product with some active constraint's normal, a constraint
 * the projection violates), so the distance is the least over the violated ones of the distance with that one added to the equalities.
 */
function nearest2(
  a: Float64Array,
  eqs: Row[],
  rows0: Row[],
  D: number,
  tolV: number,
  budget: NDBudget,
  cap: number,
): { d2: number; x: Float64Array | null } {
  // a row the whole ball of radius sqrt(cap) about a satisfies cannot matter (rows are unit vectors)
  const rad = Math.sqrt(cap);
  const rows = rows0.filter((r) => dot(r.v, 0, a, 0, D) - r.b > -rad - tolV);
  type Res = { d2: number; x: Float64Array | null };
  const NONE: Res = { d2: Infinity, x: null };
  const memo = new Map<string, Res>();
  const active: number[] = [];
  // cap: only distances below sqrt(cap) matter. Adding an equality can only lengthen the projection, and every inequality the projection
  // violates by delta (measured inside the flat) puts the nearest point at least delta farther, so a branch is dropped as soon as
  // |x0 - a|^2 + max delta^2 exceeds the cap, and a violated row is tried only while it can still beat the best found.
  const solve = (cap: number): Res => {
    const key = active
      .slice()
      .sort((x, y) => x - y)
      .join(",");
    const hit = memo.get(key);
    if (hit !== undefined) return hit;
    spend(budget);
    const E: Row[] = [...eqs, ...active.map((i) => rows[i])];
    const m = E.length;
    // Gram matrix G, its inverse by Gauss-Jordan; lambda = G^-1 r with r_i = v_i . a - b_i, x0 = a - sum lambda_i v_i
    const G: number[][] = [];
    for (let i = 0; i < m; i++) {
      G.push(new Array<number>(2 * m).fill(0));
      for (let j = 0; j < m; j++) G[i][j] = dot(E[i].v, 0, E[j].v, 0, D);
      G[i][m + i] = 1;
    }
    for (let c = 0; c < m; c++) {
      let p = c;
      for (let i = c + 1; i < m; i++) if (Math.abs(G[i][c]) > Math.abs(G[p][c])) p = i;
      // dependent equalities: this flat adds nothing
      if (Math.abs(G[p][c]) < 1e-9) {
        memo.set(key, NONE);
        return NONE;
      }
      [G[c], G[p]] = [G[p], G[c]];
      const pv = G[c][c];
      for (let j = 0; j < 2 * m; j++) G[c][j] /= pv;
      for (let i = 0; i < m; i++) {
        if (i === c) continue;
        const f = G[i][c];
        if (f !== 0) for (let j = 0; j < 2 * m; j++) G[i][j] -= f * G[c][j];
      }
    }
    const Ginv = (i: number, j: number) => G[i][m + j];
    const r: number[] = E.map((e) => dot(e.v, 0, a, 0, D) - e.b);
    const x = Float64Array.from(a);
    for (let i = 0; i < m; i++) {
      let lam = 0;
      for (let j = 0; j < m; j++) lam += Ginv(i, j) * r[j];
      for (let j = 0; j < D; j++) x[j] -= lam * E[i].v[j];
    }
    let d2 = 0;
    for (let j = 0; j < D; j++) d2 += (x[j] - a[j]) ** 2;
    if (d2 > cap) return NONE;
    const viol: [number, number][] = [];
    let worst = 0;
    const ev = new Float64Array(m);
    for (let i = 0; i < rows.length; i++) {
      if (active.includes(i)) continue;
      const v = dot(rows[i].v, 0, x, 0, D) - rows[i].b;
      if (v > tolV) {
        viol.push([v, i]);
        // delta = v / |tangential part of the row in the flat|
        for (let k = 0; k < m; k++) ev[k] = dot(E[k].v, 0, rows[i].v, 0, D);
        let q = 1;
        for (let k = 0; k < m; k++) for (let l = 0; l < m; l++) q -= ev[k] * Ginv(k, l) * ev[l];
        if (q > 1e-12) worst = Math.max(worst, (v * v) / q);
      }
    }
    if (d2 + worst > cap) return NONE;
    let out: Res = NONE;
    if (!viol.length) out = { d2, x };
    else {
      viol.sort((p, q) => q[0] - p[0]);
      for (const [, i] of viol) {
        active.push(i);
        const r2 = solve(Math.min(cap, out.d2));
        if (r2.d2 < out.d2) out = r2;
        active.pop();
      }
    }
    if (out.d2 <= cap) memo.set(key, out);
    return out;
  };
  return solve(cap);
}

const unitRow = (v: Float64Array, b: number): Row => {
  let nm = 0;
  for (let j = 0; j < v.length; j++) nm += v[j] * v[j];
  nm = Math.sqrt(nm);
  const u = new Float64Array(v.length);
  for (let j = 0; j < v.length; j++) u[j] = v[j] / nm;
  return { v: u, b: b / nm };
};

/** The rows (v . x <= b) of the window of the chain ending at node idx, in the frame of its end cell, about the source a. */
function windowRows(C: NDContainer, s: number, nodes: readonly NDNode[], idx: number): Row[] {
  const D = C.D;
  const out: Row[] = [];
  const m = new Float64Array(D);
  const aEnd = nodes[idx].a;
  // the chain's steps from the end back to the anchor; the reflections R_(end-1) ... R_p carry a vector of p's frame to the end frame
  const steps: number[] = [];
  for (let cur = idx; nodes[cur].parent >= 0; cur = nodes[cur].parent) steps.push(cur);
  // steps[0] is the end node, steps[i] the node whose parent is the (i+1)th from the end
  for (let i = 0; i < steps.length; i++) {
    const cur = steps[i];
    const p = nodes[nodes[cur].parent];
    const pc = C.cell(p.cell);
    const f = nodes[cur].via;
    const F = pc.nbr.length;
    const Sf = s * pc.hn[f] - dot(pc.nn, f * D, p.a, 0, D);
    const carry = (v: Float64Array) => {
      // apply R_p (between p and its child), then the reflections of the later steps, up to the end
      for (let t = i; t >= 0; t--) {
        const c2 = steps[t];
        const q = nodes[nodes[c2].parent];
        const mm = mirror(C.cell(q.cell), nodes[c2].via, D, m);
        reflectDir(v, m, mm, D);
      }
    };
    // heading: n_f . z >= 0  ->  -n_f . z <= 0
    {
      const v = new Float64Array(D);
      for (let j = 0; j < D; j++) v[j] = -pc.nn[f * D + j];
      carry(v);
      out.push(unitRow(v, dot(v, 0, aEnd, 0, D)));
    }
    for (let k = 0; k < F; k++) {
      if (k === f) continue;
      const Sk = s * pc.hn[k] - dot(pc.nn, k * D, p.a, 0, D);
      const v = new Float64Array(D);
      for (let j = 0; j < D; j++) v[j] = Sf * pc.nn[k * D + j] - Sk * pc.nn[f * D + j];
      carry(v);
      out.push(unitRow(v, dot(v, 0, aEnd, 0, D)));
    }
  }
  return out;
}

export interface NDHinge {
  /** Distance from the source to the nearest hinge point seen by a straight path of a chain, Infinity if none within the reach. */
  dist: number;
  /** The node and hinge pair where it occurs, and the nearest point (in that node's end-cell frame). */
  node: number;
  x: Float64Array | null;
}

/** The squared flat distance from a to the hinge flat of neighbours k, l of cell cl (within H_c): a lower bound of the distance to the hinge. */
function flatDist2(cl: NDCell, S: Float64Array, k: number, l: number, D: number): number {
  const c = dot(cl.nn, k * D, cl.nn, l * D, D) - dot(cl.nn, k * D, cl.n, 0, D) * dot(cl.nn, l * D, cl.n, 0, D);
  const det = cl.nu2[k] * cl.nu2[l] - c * c;
  if (!(det > 1e-12)) return Infinity;
  return (S[k] * S[k] * cl.nu2[l] - 2 * S[k] * S[l] * c + S[l] * S[l] * cl.nu2[k]) / det;
}

/**
 * The geodesic distance from the source to the hinges (the (d-3)-faces) over `nodes` (its half-chains), if below `reach`: for each
 * node and each hinge of its cell within reach, the nearest point of the hinge inside the chain's window.
 */
export function hingeDistance(
  C: NDContainer,
  s: number,
  nodes: readonly NDNode[],
  reach: number,
  budget: NDBudget,
  only = -1,
): NDHinge {
  const D = C.D;
  const reach2 = reach * reach * (1 + 1e-12);
  const tolV = 1e-12 * Math.max(1, s);
  let best: NDHinge = { dist: Infinity, node: -1, x: null };
  const cellEqs = (cl: NDCell): Row[] => {
    const e: Row[] = [unitRow(cl.n, s * cl.h)];
    if (C.flat) e.push(unitRow(C.flat.v, s * C.flat.b));
    return e;
  };
  for (let idx = 0; idx < nodes.length; idx++) {
    if (only >= 0 && idx !== only) continue;
    const nd = nodes[idx];
    const cl = C.cell(nd.cell);
    const F = cl.nbr.length;
    const S = new Float64Array(F);
    for (let k = 0; k < F; k++) S[k] = s * cl.hn[k] - dot(cl.nn, k * D, nd.a, 0, D);
    let wrows: Row[] | null = null;
    for (let p = 0; p < cl.hinges.length; p += 2) {
      const k = cl.hinges[p];
      const l = cl.hinges[p + 1];
      // the hinge lies in both ridge planes: farther than either plane
      if (S[k] * S[k] > reach2 * cl.nu2[k] || S[l] * S[l] > reach2 * cl.nu2[l]) continue;
      const lb = flatDist2(cl, S, k, l, D);
      if (!(lb <= reach2)) continue;
      spend(budget);
      if (!wrows) wrows = windowRows(C, s, nodes, idx);
      const eqs = [
        ...cellEqs(cl),
        unitRow(cl.nn.subarray(k * D, k * D + D), s * cl.hn[k]),
        unitRow(cl.nn.subarray(l * D, l * D + D), s * cl.hn[l]),
      ];
      const rows: Row[] = [];
      for (let j = 0; j < F; j++)
        if (j !== k && j !== l) rows.push(unitRow(cl.nn.subarray(j * D, j * D + D), s * cl.hn[j]));
      for (const r of wrows) rows.push(r);
      const r = nearest2(nd.a, eqs, rows, D, tolV, budget, Math.min(reach2, best.dist * best.dist));
      if (r.d2 < best.dist * best.dist && r.d2 <= reach2) best = { dist: Math.sqrt(r.d2), node: idx, x: r.x };
    }
  }
  return best;
}

/** The direction w (in the frame of the end cell of node idx) carried back along the chain to the anchor's frame. */
export function directionAtAnchor(
  C: NDContainer,
  nodes: readonly NDNode[],
  idx: number,
  w: Float64Array,
): Float64Array {
  const D = C.D;
  const z = Float64Array.from(w);
  const m = new Float64Array(D);
  for (let cur = idx; nodes[cur].parent >= 0; cur = nodes[cur].parent) {
    const pc = C.cell(nodes[nodes[cur].parent].cell);
    const mm = mirror(pc, nodes[cur].via, D, m);
    reflectDir(z, m, mm, D);
  }
  return z;
}
