/**
 * Platonic solids as pieces of the search (docs/spaces/BOUNDARY.md section 4C): the penalty energy and the violation of a state of rotated
 * pieces on any regular 4-polytope, and the model the search driver (search.ts) and the kernel (kernel-solids.ts) work with.
 *
 * A state holds each piece's anchor cell, its centre u (3 local coordinates, as for balls) and its orientation R (a 3 x 3 orthogonal matrix
 * in the anchor's frame, improper when a gluing reversed the handedness: every Platonic solid is mirror-symmetric, so that is the same piece
 * and `pieces` turns it into a rotation by the solid's mirror). The piece occupies { u + R p : p in the standard solid }.
 *
 * Energy: the sum of squared violations, each a penetration depth.
 *  - Pairs and loops. The chains of cells (chains.ts) are enumerated from each piece's centre within twice the solid's circumradius. A chain's
 *    development is a rigid motion T (the node's toAnchor): the piece j seen from i along that chain is the rigid copy with centre T u_j and
 *    orientation M R_j (M the linear part). When that copy's centre lies in the chain's window (the straight segment between the centres is this
 *    chain's, as for balls) the term is the squared separating-axis penetration depth of the piece i and the copy (the least overlap of their
 *    projections over the face normals and the edge-edge axes: the minimum translation that separates them). A copy of a piece itself along a
 *    non-trivial chain is a loop; it counts alike.
 *  - Edges of the polytope. For each chain within the circumradius, every hinge of the node's cell clipped to the window is a segment; the
 *    term is the squared depth of that segment inside the piece (the largest distance to the piece's faces over the segment's points).
 * This is the validator's rule (validate.ts) with one difference: the validator compares the pieces' parts clipped by the cells and windows, the
 * model the rigid copies. They agree where the development is flat, i.e. everywhere but around an edge; the final scale and nudge of the search
 * (finalizeState) is decided by the validator, so a difference can cost quality, never validity.
 *
 * Gradients are analytic in the centres and in a rotation increment of each orientation (R' = exp(theta) R, theta in the anchor's frame), by the
 * contact feature that sets the depth (a face against a vertex, or an edge against an edge); the window, the clipped segment's ends and the
 * chain's structure are held fixed. `ds`, the derivative in the edge s at fixed centres, scales the developments' translations and the hinge
 * points with s.
 */
import { ChainLimitError, clampToCell, clipSegmentToWindow, enumerateChains, inWindow } from "../chains";
import { V3, isoApply, matrixToQuaternion } from "../linalg";
import { PieceShape, pieceShape, rotationOnto } from "../pieces";
import { containerFacts, pieceBounds } from "../bounds";
import { boundaryVolume, topo } from "../polytope";
import type { BoundaryConfig, BoundaryPiece, Polytope, PolytopePieceName, Quaternion } from "../types";
import { quaternionToMatrix } from "../linalg";
import { randomCellPoint } from "./chart";
import { type PieceModel, type State, type TermOptions, copyState } from "./model";
import { randInt } from "./random";
import { M9, det9, mul9, orthonormalize, randomRotation } from "./rot";
import { clampedPieces } from "./balls";

// ---------------------------------------------------------------- the separating-axis contact and its gradient

/** The contact that sets the penetration depth of two placed solids. */
interface Contact {
  depth: number;
  n: V3;
  /** 0: a face of A, 1: a face of B, 2: an edge of each. */
  kind: number;
  eA: V3;
  eB: V3;
  L: number;
  /** The solid that sets the upper end of the overlap (0 = A, 1 = B), its vertex, and likewise for the lower end. */
  hiOwner: number;
  hiK: number;
  loOwner: number;
  loK: number;
}

const newContact = (): Contact => ({ depth: 0, n: [0, 0, 0], kind: 0, eA: [0, 0, 0], eB: [0, 0, 0], L: 1, hiOwner: 0, hiK: 0, loOwner: 0, loK: 0 });

function cross(a: V3, b: V3): V3 {
  return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
}
const dot = (a: V3, b: V3) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];

class Placed {
  readonly verts: Float64Array;
  readonly normals: Float64Array;
  readonly edges: Float64Array;
  constructor(readonly shape: PieceShape) {
    this.verts = new Float64Array(3 * shape.vertices.length);
    this.normals = new Float64Array(3 * shape.normals.length);
    this.edges = new Float64Array(3 * shape.edgeDirs.length);
  }
  set(c: V3, R: M9, ro: number): void {
    const sh = this.shape;
    for (let k = 0; k < sh.vertices.length; k++) {
      const v = sh.vertices[k];
      this.verts[3 * k] = c[0] + R[ro] * v[0] + R[ro + 1] * v[1] + R[ro + 2] * v[2];
      this.verts[3 * k + 1] = c[1] + R[ro + 3] * v[0] + R[ro + 4] * v[1] + R[ro + 5] * v[2];
      this.verts[3 * k + 2] = c[2] + R[ro + 6] * v[0] + R[ro + 7] * v[1] + R[ro + 8] * v[2];
    }
    for (let k = 0; k < sh.normals.length; k++) {
      const v = sh.normals[k];
      this.normals[3 * k] = R[ro] * v[0] + R[ro + 1] * v[1] + R[ro + 2] * v[2];
      this.normals[3 * k + 1] = R[ro + 3] * v[0] + R[ro + 4] * v[1] + R[ro + 5] * v[2];
      this.normals[3 * k + 2] = R[ro + 6] * v[0] + R[ro + 7] * v[1] + R[ro + 8] * v[2];
    }
    for (let k = 0; k < sh.edgeDirs.length; k++) {
      const v = sh.edgeDirs[k];
      this.edges[3 * k] = R[ro] * v[0] + R[ro + 1] * v[1] + R[ro + 2] * v[2];
      this.edges[3 * k + 1] = R[ro + 3] * v[0] + R[ro + 4] * v[1] + R[ro + 5] * v[2];
      this.edges[3 * k + 2] = R[ro + 6] * v[0] + R[ro + 7] * v[1] + R[ro + 8] * v[2];
    }
  }
}

/**
 * Penetration depth of two placed solids (fills `ct` with the axis that sets it). Returns the depth, <= 0 when the solids are separated
 * (touching, at 0, is allowed), stopping at the first separating axis.
 */
function satContact(A: Placed, B: Placed, ct: Contact): number {
  const nv = A.shape.vertices.length;
  let best = Infinity;
  const tryAxis = (nx: number, ny: number, nz: number, kind: number, ia: number, ib: number, L: number): boolean => {
    let maxA = -Infinity, minA = Infinity, kMaxA = 0, kMinA = 0;
    let maxB = -Infinity, minB = Infinity, kMaxB = 0, kMinB = 0;
    for (let k = 0; k < nv; k++) {
      const pa = nx * A.verts[3 * k] + ny * A.verts[3 * k + 1] + nz * A.verts[3 * k + 2];
      if (pa > maxA) { maxA = pa; kMaxA = k; }
      if (pa < minA) { minA = pa; kMinA = k; }
      const pb = nx * B.verts[3 * k] + ny * B.verts[3 * k + 1] + nz * B.verts[3 * k + 2];
      if (pb > maxB) { maxB = pb; kMaxB = k; }
      if (pb < minB) { minB = pb; kMinB = k; }
    }
    const hiOwner = maxA <= maxB ? 0 : 1;
    const loOwner = minA >= minB ? 0 : 1;
    const ov = (hiOwner === 0 ? maxA : maxB) - (loOwner === 0 ? minA : minB);
    if (ov < best) {
      best = ov;
      ct.n[0] = nx; ct.n[1] = ny; ct.n[2] = nz;
      ct.kind = kind;
      ct.hiOwner = hiOwner; ct.hiK = hiOwner === 0 ? kMaxA : kMaxB;
      ct.loOwner = loOwner; ct.loK = loOwner === 0 ? kMinA : kMinB;
      ct.L = L;
      if (kind === 2) {
        ct.eA = [A.edges[3 * ia], A.edges[3 * ia + 1], A.edges[3 * ia + 2]];
        ct.eB = [B.edges[3 * ib], B.edges[3 * ib + 1], B.edges[3 * ib + 2]];
      }
    }
    return best <= 0;
  };
  const F = A.shape.normals.length;
  for (let f = 0; f < F; f++) if (tryAxis(A.normals[3 * f], A.normals[3 * f + 1], A.normals[3 * f + 2], 0, f, 0, 1)) return (ct.depth = best);
  for (let f = 0; f < F; f++) if (tryAxis(B.normals[3 * f], B.normals[3 * f + 1], B.normals[3 * f + 2], 1, 0, f, 1)) return (ct.depth = best);
  const D = A.shape.edgeDirs.length;
  for (let a = 0; a < D; a++) {
    const ax = A.edges[3 * a], ay = A.edges[3 * a + 1], az = A.edges[3 * a + 2];
    for (let b = 0; b < D; b++) {
      const bx = B.edges[3 * b], by = B.edges[3 * b + 1], bz = B.edges[3 * b + 2];
      let nx = ay * bz - az * by, ny = az * bx - ax * bz, nz = ax * by - ay * bx;
      const L = Math.hypot(nx, ny, nz);
      if (L < 1e-7) continue;
      nx /= L; ny /= L; nz /= L;
      if (tryAxis(nx, ny, nz, 2, a, b, L)) return (ct.depth = best);
    }
  }
  return (ct.depth = best);
}

/**
 * Adds `scale` times the gradient of the contact's depth into gcA, gtA (the centre and the rotation increment of A) and gcB, gtB, by the
 * feature that sets the depth (module header). `rhoA`/`rhoB` are the vertices' offsets from their centres in the world frame.
 */
function contactGradient(A: Placed, B: Placed, cA: V3, cB: V3, ct: Contact, scale: number, gcA: V3, gtA: V3, gcB: V3, gtB: V3): void {
  const n = ct.n;
  const P = ct.hiOwner === 0 ? A : B;
  const Q = ct.loOwner === 0 ? A : B;
  const cP = ct.hiOwner === 0 ? cA : cB;
  const cQ = ct.loOwner === 0 ? cA : cB;
  const p: V3 = [P.verts[3 * ct.hiK], P.verts[3 * ct.hiK + 1], P.verts[3 * ct.hiK + 2]];
  const q: V3 = [Q.verts[3 * ct.loK], Q.verts[3 * ct.loK + 1], Q.verts[3 * ct.loK + 2]];
  const rhoP: V3 = [p[0] - cP[0], p[1] - cP[1], p[2] - cP[2]];
  const rhoQ: V3 = [q[0] - cQ[0], q[1] - cQ[1], q[2] - cQ[2]];
  const d: V3 = [p[0] - q[0], p[1] - q[1], p[2] - q[2]];
  // n . (dP - dQ): centres and the vertices' own rotation
  const gcP = ct.hiOwner === 0 ? gcA : gcB;
  const gtP = ct.hiOwner === 0 ? gtA : gtB;
  const gcQ = ct.loOwner === 0 ? gcA : gcB;
  const gtQ = ct.loOwner === 0 ? gtA : gtB;
  for (let k = 0; k < 3; k++) {
    gcP[k] += scale * n[k];
    gcQ[k] -= scale * n[k];
  }
  const rp = cross(rhoP, n);
  const rq = cross(rhoQ, n);
  for (let k = 0; k < 3; k++) {
    gtP[k] += scale * rp[k];
    gtQ[k] -= scale * rq[k];
  }
  // the axis itself moves with the solids that make it
  if (ct.kind === 0 || ct.kind === 1) {
    const nq = cross(n, d);
    const gt = ct.kind === 0 ? gtA : gtB;
    for (let k = 0; k < 3; k++) gt[k] += scale * nq[k];
  } else {
    const dn = dot(n, d);
    const dperp: V3 = [d[0] - dn * n[0], d[1] - dn * n[1], d[2] - dn * n[2]];
    const a = cross(ct.eA, cross(ct.eB, dperp));
    const b = cross(ct.eB, cross(dperp, ct.eA));
    for (let k = 0; k < 3; k++) {
      gtA[k] += (scale * a[k]) / ct.L;
      gtB[k] += (scale * b[k]) / ct.L;
    }
  }
}

// ---------------------------------------------------------------- a hinge segment inside a piece

interface HingeHit {
  depth: number;
  /** The point of the segment where the depth is largest (anchor frame). */
  x: V3;
  /** Active faces (world normals) and the weight of the first. */
  N1: V3;
  N2: V3 | null;
  w1: number;
}

/**
 * Depth of the segment a + t (b - a), t in [t0, t1], inside the placed solid (R, u): the largest, over the segment, of the least distance to
 * the faces' planes (r - max_k N_k . (x - u)). Concave piecewise linear in t: the maximum is at an end or where a rising face function crosses
 * a falling one.
 */
function hingeDepth(shape: PieceShape, normals: Float64Array, u: V3, a: V3, b: V3, t0: number, t1: number, hit: HingeHit): number {
  const F = shape.normals.length;
  const ca = new Float64Array(F);
  const cb = new Float64Array(F);
  const dx: V3 = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
  for (let k = 0; k < F; k++) {
    const nx = normals[3 * k], ny = normals[3 * k + 1], nz = normals[3 * k + 2];
    ca[k] = nx * (a[0] - u[0]) + ny * (a[1] - u[1]) + nz * (a[2] - u[2]);
    cb[k] = nx * dx[0] + ny * dx[1] + nz * dx[2];
  }
  const r = shape.inradius;
  const g = (t: number) => {
    let m = -Infinity;
    for (let k = 0; k < F; k++) m = Math.max(m, ca[k] + t * cb[k]);
    return r - m;
  };
  let bestT = t0;
  let best = g(t0);
  const g1 = g(t1);
  if (g1 > best) { best = g1; bestT = t1; }
  for (let i = 0; i < F; i++) {
    if (!(cb[i] > 0)) continue;
    for (let j = 0; j < F; j++) {
      if (!(cb[j] < 0)) continue;
      const t = (ca[j] - ca[i]) / (cb[i] - cb[j]);
      if (t >= t0 && t <= t1) {
        const v = g(t);
        if (v > best) { best = v; bestT = t; }
      }
    }
  }
  hit.depth = best;
  if (!(best > 0)) return best;
  hit.x = [a[0] + bestT * dx[0], a[1] + bestT * dx[1], a[2] + bestT * dx[2]];
  // active faces at bestT: those within a hair of the largest face function
  let m = -Infinity;
  for (let k = 0; k < F; k++) m = Math.max(m, ca[k] + bestT * cb[k]);
  const act: number[] = [];
  for (let k = 0; k < F; k++) if (ca[k] + bestT * cb[k] >= m - 1e-11) act.push(k);
  const nrm = (k: number): V3 => [normals[3 * k], normals[3 * k + 1], normals[3 * k + 2]];
  hit.N1 = nrm(act[0]);
  hit.N2 = null;
  hit.w1 = 1;
  if (bestT > t0 + 1e-12 && bestT < t1 - 1e-12 && act.length >= 2) {
    // a crossing of a rising and a falling face: the weight that cancels the t-derivative
    let rise = act[0];
    let fall = -1;
    for (const k of act) {
      if (cb[k] > cb[rise]) rise = k;
    }
    for (const k of act) if (fall < 0 || cb[k] < cb[fall]) fall = k;
    if (cb[rise] > 0 && cb[fall] < 0) {
      // lambda cb_rise + (1 - lambda) cb_fall = 0
      const lam = -cb[fall] / (cb[rise] - cb[fall]);
      hit.N1 = nrm(rise);
      hit.N2 = nrm(fall);
      hit.w1 = lam;
    }
  }
  return best;
}

// ---------------------------------------------------------------- the model

const MAX_CHAINS = 4096;

export class SolidModel implements PieceModel {
  readonly rotates = true;
  readonly target = 1e-20;
  readonly maxStep = 0.25;
  readonly floor: number;
  readonly shape: PieceShape;
  private readonly A: Placed;
  private readonly B: Placed;
  private readonly ct = newContact();

  constructor(
    readonly poly: Polytope,
    readonly piece: PolytopePieceName,
    readonly P: number,
  ) {
    this.shape = pieceShape(piece);
    this.A = new Placed(this.shape);
    this.B = new Placed(this.shape);
    this.floor = pieceBounds(poly.name, piece, P).value;
  }

  /**
   * The energy (adding its gradient into g: 3P centre entries then 3P rotation entries), or with `violation` the largest depth, over the terms
   * the options select (module header).
   */
  evaluate(cells: Int32Array, U: Float64Array, R: Float64Array, s: number, g: Float64Array | null, violation: boolean, opts?: TermOptions): number {
    const poly = this.poly;
    const T = topo(poly);
    const P = cells.length;
    const sh = this.shape;
    const W = opts?.w ?? null;
    const focus = opts?.focus ?? -1;
    const ds = opts?.ds ?? null;
    const reach = 2 * sh.circumradius;
    const byCell: number[][] = poly.cells.map(() => []);
    for (let i = 0; i < P; i++) byCell[cells[i]].push(i);
    const off = 3 * P;
    let E = 0;
    let worst = -Infinity;
    let dEds = 0;
    const gcA: V3 = [0, 0, 0], gtA: V3 = [0, 0, 0], gcB: V3 = [0, 0, 0], gtB: V3 = [0, 0, 0];
    const ct = this.ct;
    const hit: HingeHit = { depth: 0, x: [0, 0, 0], N1: [0, 0, 0], N2: null, w1: 1 };
    const Rj = new Float64Array(9);
    const nE = poly.edges.length;
    const edgeDepth = new Float64Array(nE);
    const edgeHit: (HingeHit | null)[] = new Array(nE).fill(null);
    for (let i = 0; i < P; i++) {
      if (focus >= 0 && i !== focus) continue;
      const wi = W ? W[i] : 1;
      const u = clampToCell(poly, cells[i], [U[3 * i], U[3 * i + 1], U[3 * i + 2]], s);
      let nodes;
      try {
        nodes = enumerateChains(poly, s, cells[i], u, reach, { maxChains: MAX_CHAINS });
      } catch (e) {
        if (e instanceof ChainLimitError) return Infinity;
        throw e;
      }
      this.A.set(u, R, 9 * i);
      nodes.forEach((node, ni) => {
        // ---- the edges of the polytope in this chain, inside the piece: one term per polytope edge, its deepest hit over the chains
        // (an edge on the rim of several windows is seen by each of them, and rounding decides which: the deepest does not depend on that)
        if (node.gateDist <= sh.circumradius) {
          const ct_ = T.cells[node.cell];
          const V = ct_.verts.map((p) => isoApply(node.toAnchor, [p[0] * s, p[1] * s, p[2] * s]));
          poly.cells[node.cell].hinges.forEach(([a, b], hi) => {
            const iv = clipSegmentToWindow(node.window, u, V[a], V[b]);
            if (!iv) return;
            const d = hingeDepth(sh, this.A.normals, u, V[a], V[b], iv[0], iv[1], hit);
            if (!(d > 0)) return;
            const ge = T.hingeEdge[node.cell][hi];
            if (d > edgeDepth[ge]) {
              edgeDepth[ge] = d;
              edgeHit[ge] = { depth: d, x: [hit.x[0], hit.x[1], hit.x[2]], N1: [hit.N1[0], hit.N1[1], hit.N1[2]], N2: hit.N2 ? [hit.N2[0], hit.N2[1], hit.N2[2]] : null, w1: hit.w1 };
            }
          });
        }
        // ---- the other pieces (and this one, along a non-trivial chain) in this chain's cell
        for (const j of byCell[node.cell]) {
          if ((focus < 0 && j < i) || (j === i && ni === 0)) continue;
          const uj: V3 = j === i ? u : [U[3 * j], U[3 * j + 1], U[3 * j + 2]];
          const y = isoApply(node.toAnchor, uj);
          if (!inWindow(node, u, y)) continue;
          const dist = Math.hypot(y[0] - u[0], y[1] - u[1], y[2] - u[2]);
          if (dist >= reach) continue;
          if (j === i && dist <= 1e-12) continue;
          mul9(node.toAnchor.m, 0, R, 9 * j, Rj, 0);
          this.B.set(y, Rj, 0);
          const depth = satContact(this.A, this.B, ct);
          if (!(depth > 0)) continue;
          if (violation) {
            worst = Math.max(worst, depth);
            continue;
          }
          const wij = W ? 0.5 * (W[i] + W[j]) : 1;
          E += wij * depth * depth;
          if (!g && !ds) continue;
          gcA[0] = gcA[1] = gcA[2] = gtA[0] = gtA[1] = gtA[2] = 0;
          gcB[0] = gcB[1] = gcB[2] = gtB[0] = gtB[1] = gtB[2] = 0;
          contactGradient(this.A, this.B, u, y, ct, 2 * wij * depth, gcA, gtA, gcB, gtB);
          const m = node.toAnchor.m;
          const detm = det9(m);
          if (ds) dEds += (gcB[0] * node.toAnchor.t[0] + gcB[1] * node.toAnchor.t[1] + gcB[2] * node.toAnchor.t[2]) / s;
          if (!g) continue;
          for (let k = 0; k < 3; k++) {
            g[3 * i + k] += gcA[k];
            g[off + 3 * i + k] += gtA[k];
            // back to piece j's own frame: M^T for the centre, det(M) M^T for the rotation increment
            g[3 * j + k] += m[k] * gcB[0] + m[3 + k] * gcB[1] + m[6 + k] * gcB[2];
            g[off + 3 * j + k] += detm * (m[k] * gtB[0] + m[3 + k] * gtB[1] + m[6 + k] * gtB[2]);
          }
        }
      });
      // the edges' terms
      for (let e = 0; e < nE; e++) {
        const d = edgeDepth[e];
        if (!(d > 0)) continue;
        edgeDepth[e] = 0;
        if (violation) {
          worst = Math.max(worst, d);
          continue;
        }
        E += wi * d * d;
        if (!g && !ds) continue;
        const h = edgeHit[e]!;
        const c = 2 * wi * d;
        const dxu: V3 = [h.x[0] - u[0], h.x[1] - u[1], h.x[2] - u[2]];
        // depth = r - N . (x - u) at the active face(s): d/du = N, d/dtheta = (x - u) x N
        const faces: [V3, number][] = h.N2 ? [[h.N1, h.w1], [h.N2, 1 - h.w1]] : [[h.N1, 1]];
        for (const [N, w] of faces) {
          if (g) {
            const gt = cross(dxu, N);
            for (let k = 0; k < 3; k++) {
              g[3 * i + k] += c * w * N[k];
              g[off + 3 * i + k] += c * w * gt[k];
            }
          }
          if (ds) dEds -= (c * w * dot(N, h.x)) / s;
        }
      }
    }
    if (ds) ds[0] = dEds;
    return violation ? worst : E;
  }

  energy(st: State, s: number, g: Float64Array): number {
    g.fill(0);
    return this.evaluate(st.cells, st.U, st.R!, s, g, false);
  }

  violation(st: State, s: number): number {
    return this.evaluate(st.cells, st.U, st.R!, s, null, true);
  }

  terms(st: State, s: number, g: Float64Array | null, violation: boolean, o: TermOptions): number {
    if (g) g.fill(0);
    if (o.ds) o.ds[0] = 0;
    return this.evaluate(st.cells, st.U, st.R!, s, g, violation, o);
  }

  /**
   * A random start edge: pieces at a quarter of the boundary's volume, never below 1.35 of the proven bound, and roomy enough that a cell's
   * inscribed ball holds the solid's circumscribed ball (a start whose pieces each span many cells and edges is hard to relax).
   */
  startEdge(): number {
    const f = this.floor;
    const vol = Math.cbrt((4 * this.P * this.shape.volume) / boundaryVolume(this.poly));
    const room = this.shape.circumradius / containerFacts(this.poly.name).inradius;
    return Math.max(1.35 * f, vol, room);
  }

  randomState(s: number, rand: () => number): State {
    const C = this.poly.cells.length;
    const cells = new Int32Array(this.P);
    const U = new Float64Array(3 * this.P);
    const R = new Float64Array(9 * this.P);
    for (let i = 0; i < this.P; i++) {
      cells[i] = randInt(rand, C);
      U.set(randomCellPoint(this.poly, cells[i], s, rand, 0.6), 3 * i);
      R.set(randomRotation(rand), 9 * i);
    }
    return { cells, U, R };
  }

  hop(st: State, s: number, rand: () => number): State {
    const out = copyState(st);
    const k = 1 + randInt(rand, 3);
    for (let t = 0; t < k; t++) {
      const i = randInt(rand, this.P);
      out.cells[i] = randInt(rand, this.poly.cells.length);
      out.U.set(randomCellPoint(this.poly, out.cells[i], s, rand, 0.6), 3 * i);
      out.R!.set(randomRotation(rand), 9 * i);
    }
    return out;
  }

  /**
   * The cells themselves, when the solid is the polytope's own cell and P <= C: each piece centred in a cell and turned onto it, valid from
   * s = 1 (a hair above); null otherwise.
   */
  constructedStart(rand: () => number): { st: State; s: number } | null {
    const C = this.poly.cells.length;
    if (this.P > C) return null;
    const Rs: number[][] = [];
    for (const c of this.poly.cells) {
      const R = rotationOnto(this.shape, c.vertices as V3[]);
      if (!R) return null;
      Rs.push(R);
    }
    const order = Array.from({ length: C }, (_, i) => i);
    const cells = new Int32Array(this.P);
    const R = new Float64Array(9 * this.P);
    for (let i = 0; i < this.P; i++) {
      const j = i + randInt(rand, C - i);
      [order[i], order[j]] = [order[j], order[i]];
      cells[i] = order[i];
      R.set(Rs[order[i]], 9 * i);
    }
    return { st: { cells, U: new Float64Array(3 * this.P), R }, s: 1 + 1e-6 };
  }

  stateOf(config: BoundaryConfig): State | null {
    if (!config || !Array.isArray(config.pieces) || config.pieces.length !== this.P) return null;
    if (!(typeof config.s === "number" && Number.isFinite(config.s) && config.s > 0)) return null;
    const cells = new Int32Array(this.P);
    const U = new Float64Array(3 * this.P);
    const R = new Float64Array(9 * this.P);
    for (let i = 0; i < this.P; i++) {
      const p = config.pieces[i];
      if (!p || !Number.isInteger(p.cell) || p.cell < 0 || p.cell >= this.poly.cells.length || !Array.isArray(p.u) || p.u.length !== 3 || !p.u.every(Number.isFinite)) return null;
      cells[i] = p.cell;
      U.set(clampToCell(this.poly, p.cell, p.u, config.s), 3 * i);
      R.set(p.q ? quaternionToMatrix(p.q) : [1, 0, 0, 0, 1, 0, 0, 0, 1], 9 * i);
    }
    return { cells, U, R };
  }

  pieces(st: State, s: number): BoundaryPiece[] {
    const out = clampedPieces(this.poly, st, s);
    const S = this.shape.mirror;
    const tmp = new Float64Array(9);
    for (let i = 0; i < out.length; i++) {
      const Ri = Array.from(st.R!.subarray(9 * i, 9 * i + 9));
      orthonormalize(Ri);
      // an improper orientation is the same piece turned by the solid's mirror
      if (det9(Ri) < 0) {
        mul9(Ri, 0, S, 0, tmp, 0);
        for (let k = 0; k < 9; k++) Ri[k] = tmp[k];
        orthonormalize(Ri);
      }
      out[i].q = matrixToQuaternion(Ri) as Quaternion;
    }
    return out;
  }
}
