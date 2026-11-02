/**
 * The penalty model of balls on the simplex, orthoplex and demicube (d = 5..10), the PieceModel of search/model.ts for the nD
 * containers; its engine is chains.ts (docs/spaces/BOUNDARY.md 4D). A state holds each ball's anchor cell and centre (n = d - 1 local
 * coordinates at the state's edge).
 *
 * Energy = sum over balls of w_i (1/2 - h_i)^2, h_i the geodesic distance of the centre to the hinges when it is below 1/2, plus sum
 * over pairs of w_ij (1 - l_ij)^2, l_ij the length of the shortest straight path joining the centres when it is below 1, plus the
 * same for each geodesic loop at a centre (w_ij = (w_i + w_j) / 2). The gradient is analytic in the centres: a path is a straight
 * segment between a centre and the image a = T^(-1) x of the other, so d l / d x_i = L (a_i - a_j) / l with L the linear part of the
 * development; the hinge distance h = |x* - a| has d h / d x = L (a - x*) / h at the nearest point x* (exact when no window plane is
 * active at x*, a descent direction otherwise: the finishing stage judges by the exact rule). The edge-derivative at fixed centres is
 * a central difference (the kernel adds the centres' own motion). A ball nearer than 1/2 to a hinge of its own cell is not developed:
 * its pair terms fall back to the chord in R^D (a lower bound of the geodesic), so a crowded start costs little.
 */
import type { BoundaryConfig, BoundaryPiece } from "../types";
import type { PieceModel, State, TermOptions } from "../search/model";
import { ndBound, ndCentreEdge, ndPairEdge } from "./bounds";
import { type NDNode, NDLimitError, develop, hingeDistance, joinDistance, ndZeroLength, newNDBudget } from "./chains";
import { type NDContainer, getND } from "./container";
import { type NDName, ndParse } from "./names";
import { type NDAnchored, ND_REACH, anchorPiece, ndSnap } from "./validate";

/** Steps one evaluation of the model may spend before it reports a huge violation instead. */
export const ND_EVAL_WORK = 60_000;
const HUGE = 1e6;

export class NDBallModel implements PieceModel {
  readonly target = 1e-22;
  readonly maxStep = 0.5;
  readonly floor: number;
  readonly n: number;
  readonly C: NDContainer;

  constructor(
    readonly polytope: string,
    readonly P: number,
  ) {
    const p = ndParse(polytope);
    if (!p) throw new Error(`NDBallModel: ${polytope} is none of the simplex, orthoplex or demicube containers`);
    this.C = getND(polytope as Parameters<typeof getND>[0]);
    this.n = this.C.n;
    this.floor = ndBound(p.family, p.d, P);
  }

  /** Evaluates the terms; g (when given) receives the gradient in U, `viol` makes it return the largest violation instead of the energy. */
  private evaluate(st: State, s: number, g: Float64Array | null, viol: boolean, opts?: TermOptions): number {
    const { C, P, n } = this;
    const D = C.D;
    const W = opts?.w ?? null;
    const focus = opts?.focus ?? -1;
    const wantGrad = g !== null;
    const budget = newNDBudget(ND_EVAL_WORK);
    const eta = ndSnap(s, 0);
    let E = 0;
    let worst = -Infinity;
    if (g) g.fill(0);
    const amb = Array.from({ length: P }, () => (wantGrad ? new Float64Array(D) : null));
    try {
      const at: (NDAnchored | undefined)[] = new Array(P);
      const chains: (NDNode[][] | undefined)[] = new Array(P);
      const hinge: number[] = new Array(P).fill(Infinity);
      const crowded = new Uint8Array(P);
      const need = new Uint8Array(P);
      const pos = (i: number) => (at[i] ??= anchorPiece(C, s, st.cells[i], st.U, eta, i * n));
      // chord-near pairs first, to know which balls the focus needs
      const near: [number, number][] = [];
      for (let i = 0; i < P; i++) {
        if (focus >= 0 && i !== focus) continue;
        const xi = pos(i).x;
        need[i] = 1;
        for (let j = 0; j < P; j++) {
          if (j === i || (focus < 0 && j < i)) continue;
          const xj = pos(j).x;
          let q = 0;
          for (let r = 0; r < D; r++) q += (xi[r] - xj[r]) ** 2;
          if (q < 1) {
            near.push([i, j]);
            need[j] = 1;
          }
        }
      }
      // clearance, and the half-chains of the balls that are not crowded against a hinge
      for (let i = 0; i < P; i++) {
        if (!need[i]) continue;
        const a = pos(i);
        const roots = a.anchors.map((an) => develop(C, s, an.cell, an.x, 0, budget));
        let best = { dist: Infinity, node: -1, x: null as Float64Array | null, set: -1 };
        for (let k = 0; k < roots.length; k++) {
          const h = hingeDistance(C, s, roots[k], ND_REACH, budget);
          if (h.dist < best.dist) best = { dist: h.dist, node: h.node, x: h.x, set: k };
        }
        let chainsOf: NDNode[][] | undefined;
        if (best.dist >= 0.5) {
          chainsOf = a.anchors.map((an) => develop(C, s, an.cell, an.x, ND_REACH, budget));
          chains[i] = chainsOf;
          for (let k = 0; k < chainsOf.length; k++) {
            const h = hingeDistance(C, s, chainsOf[k], ND_REACH, budget);
            if (h.dist < best.dist) best = { dist: h.dist, node: h.node, x: h.x, set: k };
          }
        } else crowded[i] = 1;
        hinge[i] = best.dist;
        const wi = W ? W[i] : 1;
        if (focus >= 0 && i !== focus) continue;
        const v = 0.5 - best.dist;
        if (viol) worst = Math.max(worst, v);
        else if (v > 0) {
          E += wi * v * v;
          if (amb[i] && best.x && best.dist > 0) {
            const nodeSet = crowded[i] ? roots[best.set] : chainsOf![best.set];
            const nd = nodeSet[best.node];
            const c = (-2 * wi * v) / best.dist;
            // d h / d x = L (a - x*) / h
            const diff = new Float64Array(D);
            for (let r = 0; r < D; r++) diff[r] = nd.a[r] - best.x[r];
            for (let r = 0; r < D; r++) {
              let t = 0;
              for (let q = 0; q < D; q++) t += nd.L[r * D + q] * diff[q];
              amb[i]![r] += c * t;
            }
          }
        }
      }
      // pairs
      for (const [i, j] of near) {
        const xi = at[i]!.x;
        const xj = at[j]!.x;
        const wij = W ? 0.5 * (W[i] + W[j]) : 1;
        let ell = Infinity;
        let ga: Float64Array | null = null;
        let gb: Float64Array | null = null;
        let q0 = 0;
        for (let r = 0; r < D; r++) q0 += (xi[r] - xj[r]) ** 2;
        if (crowded[i] || crowded[j] || q0 < ndZeroLength(s) ** 2) {
          let q = 0;
          for (let r = 0; r < D; r++) q += (xi[r] - xj[r]) ** 2;
          ell = Math.sqrt(q);
          if (wantGrad && ell > 0) {
            ga = new Float64Array(D);
            gb = new Float64Array(D);
            for (let r = 0; r < D; r++) {
              ga[r] = (xi[r] - xj[r]) / ell;
              gb[r] = -ga[r];
            }
          }
        } else {
          let bestJ: { dist: number; a: number; b: number } = { dist: Infinity, a: -1, b: -1 };
          let bi = -1;
          let bj = -1;
          for (let p = 0; p < chains[i]!.length; p++)
            for (let q = 0; q < chains[j]!.length; q++) {
              const r = joinDistance(
                C,
                s,
                chains[i]![p],
                chains[j]![q],
                Math.min(bestJ.dist, 1),
                budget,
                false,
                1e-12,
                false,
              );
              if (r.dist < bestJ.dist) {
                bestJ = r;
                bi = p;
                bj = q;
              }
            }
          ell = bestJ.dist;
          if (wantGrad && ell < 1 && ell > 0) {
            const A = chains[i]![bi][bestJ.a];
            const B = chains[j]![bj][bestJ.b];
            ga = new Float64Array(D);
            gb = new Float64Array(D);
            for (let r = 0; r < D; r++) {
              let ta = 0;
              let tb = 0;
              for (let q = 0; q < D; q++) {
                ta += A.L[r * D + q] * (A.a[q] - B.a[q]);
                tb += B.L[r * D + q] * (B.a[q] - A.a[q]);
              }
              ga[r] = ta / ell;
              gb[r] = tb / ell;
            }
          }
        }
        if (!(ell < 1)) continue;
        const v = 1 - ell;
        if (viol) worst = Math.max(worst, v);
        else {
          E += wij * v * v;
          if (ga && gb) {
            const c = -2 * wij * v;
            for (let r = 0; r < D; r++) {
              if (amb[i]) amb[i]![r] += c * ga[r];
              if (amb[j]) amb[j]![r] += c * gb[r];
            }
          }
        }
      }
      // loops
      for (let i = 0; i < P; i++) {
        if (!need[i] || (focus >= 0 && i !== focus) || crowded[i] || !chains[i]) continue;
        const wi = W ? W[i] : 1;
        let bestL = { dist: Infinity, a: -1, b: -1 };
        let bp = -1;
        let bq = -1;
        for (let p = 0; p < chains[i]!.length; p++)
          for (let q = p; q < chains[i]!.length; q++) {
            const r = joinDistance(C, s, chains[i]![p], chains[i]![q], Math.min(bestL.dist, 1), budget, p === q);
            if (r.dist < bestL.dist) {
              bestL = r;
              bp = p;
              bq = q;
            }
          }
        if (!(bestL.dist < 1)) continue;
        const v = 1 - bestL.dist;
        if (viol) worst = Math.max(worst, v);
        else {
          E += wi * v * v;
          if (amb[i] && bestL.dist > 0) {
            const A = chains[i]![bp][bestL.a];
            const B = chains[i]![bq][bestL.b];
            const c = (-2 * wi * v) / bestL.dist;
            for (let r = 0; r < D; r++) {
              let t = 0;
              for (let q = 0; q < D; q++) t += A.L[r * D + q] * (A.a[q] - B.a[q]) + B.L[r * D + q] * (B.a[q] - A.a[q]);
              amb[i]![r] += c * t;
            }
          }
        }
      }
    } catch (e) {
      if (!(e instanceof NDLimitError)) throw e;
      return viol ? 1 : HUGE;
    }
    if (viol) return worst === -Infinity ? -1 : worst;
    if (g) {
      // ambient gradient -> local coordinates of the state's cell
      for (let i = 0; i < P; i++) {
        const a = amb[i];
        if (!a) continue;
        const fr = C.cell(st.cells[i]).frame;
        for (let k = 0; k < n; k++) {
          let t = 0;
          for (let r = 0; r < D; r++) t += fr[k][r] * a[r];
          g[i * n + k] += t;
        }
      }
      if (opts?.ds) {
        const h = 1e-6 * s;
        const up = this.evaluate(st, s + h, null, false, { w: opts.w, focus: opts.focus });
        const dn = this.evaluate(st, s - h, null, false, { w: opts.w, focus: opts.focus });
        opts.ds[0] = (up - dn) / (2 * h);
      }
    }
    return E;
  }

  energy(st: State, s: number, g: Float64Array): number {
    return this.evaluate(st, s, g, false);
  }

  violation(st: State, s: number): number {
    return this.evaluate(st, s, null, true);
  }

  terms(st: State, s: number, g: Float64Array | null, violation: boolean, o: TermOptions): number {
    return this.evaluate(st, s, g, violation, o);
  }

  startEdge(): number {
    // a ball fits at a cell's centre from ndCentreEdge, which is above the bound on the demicube (d >= 6); start 30% above it
    return Math.max(
      this.floor * 1.5,
      ndCentreEdge(this.polytope as NDName) * 1.3,
      this.P > 1 ? ndPairEdge(this.polytope as NDName) * 1.1 : 0,
    );
  }

  /** A random point of cell c (local coordinates at edge s), shrunk toward the cell's centre until it is inside. */
  randomPoint(c: number, s: number, rand: () => number, spread = 0.2 / Math.sqrt(this.n)): number[] {
    const { C, n } = this;
    const cl = C.cell(c);
    const D = C.D;
    const u0 = new Array<number>(n);
    for (let k = 0; k < n; k++)
      u0[k] = Math.sqrt(-2 * Math.log(1 - rand())) * Math.cos(2 * Math.PI * rand()) * spread * s;
    let lam = 1;
    for (let t = 0; t < 40; t++) {
      const u = u0.map((v) => v * lam);
      const x = C.toAmbient(c, u, s);
      let ok = true;
      for (let k = 0; k < cl.nbr.length && ok; k++) {
        let v = -s * cl.hn[k];
        for (let r = 0; r < D; r++) v += cl.nn[k * D + r] * x[r];
        if (v > 0) ok = false;
      }
      if (ok) return u;
      lam *= 0.7;
    }
    return new Array<number>(n).fill(0);
  }

  /** A random cell a ball can sit in at edge s: the demicube's simplex facets (clearance s / sqrt(n^2 - 1) at the centre) only when large enough. */
  viableCell(rand: () => number, s: number): number {
    const { C, n } = this;
    for (let t = 0; t < 30; t++) {
      const c = C.pickCell(rand);
      if (C.family !== "demicube" || c < 2 * C.d || s >= 0.55 * Math.sqrt(n * n - 1)) return c;
    }
    return 0;
  }

  randomState(s: number, rand: () => number): State {
    const { C, P, n } = this;
    const cells = new Int32Array(P);
    const U = new Float64Array(n * P);
    // distinct cells while there are enough viable ones (a ball per cell is how these containers hold the fewest balls)
    const pool: number[] = [];
    const nPool = C.family === "demicube" && s < 0.55 * Math.sqrt(n * n - 1) ? 2 * C.d : C.C;
    if (nPool <= 4096) for (let c = 0; c < nPool; c++) pool.push(c);
    for (let i = 0; i < P; i++) {
      if (i < pool.length) {
        const j = i + Math.min(pool.length - i - 1, Math.floor(rand() * (pool.length - i)));
        [pool[i], pool[j]] = [pool[j], pool[i]];
        cells[i] = pool[i];
      } else cells[i] = this.viableCell(rand, s);
      U.set(this.randomPoint(cells[i], s, rand), i * n);
    }
    return { cells, U };
  }

  hop(st: State, s: number, rand: () => number): State {
    const { C, P, n } = this;
    const out: State = { cells: Int32Array.from(st.cells), U: Float64Array.from(st.U) };
    const m = 1 + Math.floor(rand() * 3);
    for (let t = 0; t < m; t++) {
      const i = Math.floor(rand() * P);
      out.cells[i] = this.viableCell(rand, s);
      out.U.set(this.randomPoint(out.cells[i], s, rand), i * n);
    }
    return out;
  }

  constructedStart(): { st: State; s: number } | null {
    return null;
  }

  stateOf(config: BoundaryConfig): State | null {
    const { C, P, n } = this;
    if (!config || !Array.isArray(config.pieces) || config.pieces.length !== P) return null;
    if (!(typeof config.s === "number" && Number.isFinite(config.s) && config.s > 0)) return null;
    const cells = new Int32Array(P);
    const U = new Float64Array(n * P);
    for (let i = 0; i < P; i++) {
      const p = config.pieces[i];
      if (
        !p ||
        !Number.isInteger(p.cell) ||
        p.cell < 0 ||
        p.cell >= C.C ||
        !Array.isArray(p.u) ||
        p.u.length !== n ||
        !p.u.every(Number.isFinite)
      )
        return null;
      cells[i] = p.cell;
      for (let k = 0; k < n; k++) U[i * n + k] = p.u[k];
    }
    return { cells, U };
  }

  /** The pieces of a state: each centre pulled toward its cell's centre by a hair until its reading is inside the closed cell. */
  pieces(st: State, s: number): BoundaryPiece[] {
    const { C, P, n } = this;
    const D = C.D;
    const out: BoundaryPiece[] = [];
    for (let i = 0; i < P; i++) {
      const c = st.cells[i];
      const cl = C.cell(c);
      let u = Array.from(st.U.subarray(i * n, i * n + n));
      for (let t = 0, f = 1e-16; t < 60; t++, f *= 2) {
        const x = C.toAmbient(c, u, s);
        let worst = -Infinity;
        for (let k = 0; k < cl.nbr.length; k++) {
          let v = -s * cl.hn[k];
          for (let r = 0; r < D; r++) v += cl.nn[k * D + r] * x[r];
          if (v > worst) worst = v;
        }
        if (worst <= 0) break;
        u = Array.from(st.U.subarray(i * n, i * n + n), (v) => v * (1 - f));
      }
      out.push({ cell: c, u });
    }
    return out;
  }
}
