/**
 * The penalty models of the d-cube kernel (kernel.ts): balls (diameter 1) and axis-aligned unit cubes, the d-cube twins of
 * search/balls.ts (WordMetric) and search/cubes.ts (AlignedCubeModel), written for n = d - 1 local axes. A state holds each
 * piece's anchor cell and centre (n numbers, at the state's edge s); the models are the PieceModel of search/model.ts, so
 * the search's finalize stage (finalizeState) serves them.
 *
 * Balls. The energy is the sum of the squared violations: (1/2 - clearance)^2 for every pair of axes whose nearest hinge is
 * nearer than 1/2 (the clearance has the closed form of lattice.ts, one term per pair of axes), and (1 - d)^2 for every
 * realised straight geodesic of length d < 1 between two centres, found by reach (lattice.ts realisedWords). A geodesic
 * loop at a centre is not a term: a word that ends in its start cell crosses two lattice axes twice (BOUNDARY.md 4B), so a
 * loop is at least s sqrt 2 long, which is >= 1 once a ball fits (s >= 1/sqrt 2). The gradient is analytic: a geodesic is
 * a straight segment from u_i to the image y = M u_j + h t, so d(d)/du_i = (u_i - y)/d and d(d)/du_j = M^T (y - u_i)/d.
 *
 * Cubes (edge 1, axis aligned). The cube centred at u is the box u +- 1/2 and crosses axis k when |u_k| + 1/2 > h; its parts
 * are the box clipped to the anchor and, for each crossed axis, the slab beyond the face, in the neighbour. The energy is the
 * squared crossing depths except the largest (a hinge inside the cube), plus, for every two parts of distinct cubes in a
 * cell, the squared penetration depth (the least overlap over the axes). Every box bound is a u + b in a single coordinate,
 * so the gradient is exact.
 */
import type { BoundaryConfig, BoundaryPiece } from "../types";
import type { PieceModel, State, TermOptions } from "../search/model";
import { randInt } from "../search/random";
import { ncubeBound } from "./bounds";
import { nearPairs } from "./validate";
import { WordLimitError, axesOf, embedCell, newBudget, realisedWords, type Realised } from "./lattice";

/** A random point of the cell (n coordinates), `shrink` times the way from its centre to its faces. */
export function randomPoint(n: number, s: number, rand: () => number, shrink: number): number[] {
  const h = (s / 2) * shrink;
  return Array.from({ length: n }, () => (2 * rand() - 1) * h);
}

/** Cells and centres of a configuration with P pieces in range, or null; centres are clamped into their closed cells. */
export function plainState(d: number, config: BoundaryConfig, P: number): State | null {
  const n = d - 1;
  if (!config || !Array.isArray(config.pieces) || config.pieces.length !== P) return null;
  if (!(typeof config.s === "number" && Number.isFinite(config.s) && config.s > 0)) return null;
  const h = config.s / 2;
  const cells = new Int32Array(P);
  const U = new Float64Array(n * P);
  for (let i = 0; i < P; i++) {
    const p = config.pieces[i];
    if (
      !p ||
      !Number.isInteger(p.cell) ||
      p.cell < 0 ||
      p.cell >= 2 * d ||
      !Array.isArray(p.u) ||
      p.u.length !== n ||
      !p.u.every(Number.isFinite)
    )
      return null;
    cells[i] = p.cell;
    for (let k = 0; k < n; k++) U[i * n + k] = Math.max(-h, Math.min(h, p.u[k]));
  }
  return { cells, U };
}

/** The pieces of a state, each centre clamped into its closed cell exactly (the validator checks it at tol 0). */
export function clampedPieces(d: number, st: State, s: number): BoundaryPiece[] {
  const n = d - 1;
  const h = s * 0.5;
  const out: BoundaryPiece[] = [];
  for (let i = 0; i < st.cells.length; i++) {
    const u: number[] = [];
    for (let k = 0; k < n; k++) u.push(Math.min(h, Math.max(-h, st.U[i * n + k])));
    out.push({ cell: st.cells[i], u });
  }
  return out;
}

// ---------------------------------------------------------------- balls

/** Search steps one pair of balls may use in the energy before the chord stands in for its geodesics (a valid layout needs a few hundred). */
export const PAIR_WORK = 8_000;

export class NCubeBallModel implements PieceModel {
  readonly target = 1e-22;
  readonly maxStep = 0.5;
  readonly floor: number;
  readonly n: number;

  constructor(
    readonly d: number,
    readonly P: number,
  ) {
    this.n = d - 1;
    this.floor = ncubeBound(d, "ball", P);
  }

  private evaluate(st: State, s: number, g: Float64Array | null, violation: boolean, opts?: TermOptions): number {
    const { d, n, P } = this;
    const h = s / 2;
    const W = opts?.w ?? null;
    const focus = opts?.focus ?? -1;
    const ds = opts?.ds ?? null;
    const U = st.U;
    const cells = st.cells;
    let dEds = 0;
    let E = 0;
    let worst = -Infinity;
    // clearance, closed form: one term per pair of axes
    for (let i = 0; i < P; i++) {
      if (focus >= 0 && i !== focus) continue;
      const wi = W ? W[i] : 1;
      const o = n * i;
      for (let a = 0; a < n; a++) {
        const Aa = h - Math.abs(U[o + a]);
        if (Aa >= 0.5 && !violation) continue;
        for (let b = a + 1; b < n; b++) {
          const Ab = h - Math.abs(U[o + b]);
          const cl = Math.sqrt(Aa * Aa + Ab * Ab);
          if (violation) {
            if (0.5 - cl > worst) worst = 0.5 - cl;
          } else if (cl < 0.5) {
            E += wi * (0.5 - cl) * (0.5 - cl);
            if ((g || ds) && cl > 0) {
              const c = (2 * wi * (0.5 - cl)) / cl;
              if (g) {
                g[o + a] += c * Aa * (U[o + a] >= 0 ? 1 : -1);
                g[o + b] += c * Ab * (U[o + b] >= 0 ? 1 : -1);
              }
              if (ds) dEds -= 0.5 * c * (Aa + Ab);
            }
          }
        }
      }
    }
    // pairs whose R^d centres are < 1 apart (a geodesic is never shorter than the chord)
    const pairs =
      P <= 48
        ? allPairs(d, cells, U, n, h, P, focus)
        : nearPairs(d, Array.from(cells), U, n, h, 1).filter(([i, j]) => focus < 0 || i === focus || j === focus);
    for (const [i, j] of pairs) {
      const wij = W ? 0.5 * (W[i] + W[j]) : 1;
      let found: Realised[];
      try {
        found = realisedWords(d, s, cells[i], U, i * n, cells[j], U, j * n, newBudget(PAIR_WORK));
      } catch (e) {
        if (!(e instanceof WordLimitError)) throw e;
        // a pair buried in each other near the hinges of a container far too small for them (s well below 1/sqrt 2): the
        // search for its geodesics is cut off, and the chord in R^d, a lower bound of every geodesic, stands in for it
        const r = this.chordTerm(cells[i], cells[j], U, i, j, h, wij, g, violation, ds);
        if (violation) worst = Math.max(worst, r.violation);
        else {
          E += r.energy;
          dEds += r.dEds;
        }
        continue;
      }
      for (const r of found) {
        const dist = r.len;
        if (violation) {
          if (1 - dist > worst) worst = 1 - dist;
          continue;
        }
        if (dist >= 1) continue;
        E += wij * (1 - dist) * (1 - dist);
        if (!g && !ds) continue;
        const c = (2 * wij * (1 - dist)) / Math.max(dist, 1e-15);
        for (let k = 0; k < n; k++) {
          const rk = r.y[k] - U[i * n + k];
          if (ds) dEds -= 0.5 * c * rk * r.t[k];
          if (!g) continue;
          // dE/du_i = c r and dE/du_j = -c M^T r, r = y - u_i (as search/balls.ts)
          g[i * n + k] += c * rk;
          g[j * n + r.perm[k]] -= c * r.sign[k] * rk;
        }
      }
    }
    if (ds) ds[0] += dEds;
    return violation ? worst : E;
  }

  /** The penalty of a pair by the chord in R^d (a lower bound of the geodesic), with its gradient and dE/ds. */
  private chordTerm(
    ci: number,
    cj: number,
    U: Float64Array,
    i: number,
    j: number,
    h: number,
    w: number,
    g: Float64Array | null,
    violation: boolean,
    ds: Float64Array | null,
  ): { energy: number; dEds: number; violation: number } {
    const { d, n } = this;
    const xi = embedCell(d, ci, U, i * n, h);
    const xj = embedCell(d, cj, U, j * n, h);
    let q = 0;
    for (let t = 0; t < d; t++) q += (xi[t] - xj[t]) ** 2;
    const chord = Math.sqrt(q);
    if (violation) return { energy: 0, dEds: 0, violation: 1 - chord };
    if (!(chord < 1)) return { energy: 0, dEds: 0, violation: -Infinity };
    const c = (2 * w * (1 - chord)) / Math.max(chord, 1e-15);
    let dEds = 0;
    const axi = axesOf(d, ci >> 1);
    const axj = axesOf(d, cj >> 1);
    if (g) {
      for (let k = 0; k < n; k++) {
        g[i * n + k] -= c * (xi[axi[k]] - xj[axi[k]]);
        g[j * n + k] += c * (xi[axj[k]] - xj[axj[k]]);
      }
    }
    if (ds) {
      // only the coordinate along each normal moves with s (sigma h)
      dEds -=
        c * ((xi[ci >> 1] - xj[ci >> 1]) * (ci & 1 ? 0.5 : -0.5) - (xi[cj >> 1] - xj[cj >> 1]) * (cj & 1 ? 0.5 : -0.5));
    }
    return { energy: w * (1 - chord) * (1 - chord), dEds, violation: -Infinity };
  }

  energy(st: State, s: number, g: Float64Array): number {
    g.fill(0);
    return this.evaluate(st, s, g, false);
  }

  violation(st: State, s: number): number {
    return this.evaluate(st, s, null, true);
  }

  terms(st: State, s: number, g: Float64Array | null, violation: boolean, o: TermOptions): number {
    if (g) g.fill(0);
    if (o.ds) o.ds[0] = 0;
    return this.evaluate(st, s, g, violation, o);
  }

  startEdge(): number {
    return 1.35 * this.floor;
  }

  randomState(s: number, rand: () => number): State {
    const { d, n, P } = this;
    const cells = new Int32Array(P);
    const U = new Float64Array(n * P);
    for (let i = 0; i < P; i++) {
      cells[i] = randInt(rand, 2 * d);
      U.set(randomPoint(n, s, rand, 0.6), n * i);
    }
    return { cells, U };
  }

  hop(st: State, s: number, rand: () => number): State {
    const { d, n, P } = this;
    const out: State = { cells: Int32Array.from(st.cells), U: Float64Array.from(st.U) };
    const k = 1 + randInt(rand, 3);
    for (let t = 0; t < k; t++) {
      const i = randInt(rand, P);
      out.cells[i] = randInt(rand, 2 * d);
      out.U.set(randomPoint(n, s, rand, 0.6), n * i);
    }
    return out;
  }

  constructedStart(): null {
    return null;
  }

  stateOf(config: BoundaryConfig): State | null {
    return plainState(this.d, config, this.P);
  }

  pieces(st: State, s: number): BoundaryPiece[] {
    return clampedPieces(this.d, st, s);
  }
}

function allPairs(
  d: number,
  cells: Int32Array,
  U: Float64Array,
  n: number,
  h: number,
  P: number,
  focus: number,
): [number, number][] {
  const X: Float64Array[] = [];
  for (let i = 0; i < P; i++) X.push(embedCell(d, cells[i], U, i * n, h));
  const out: [number, number][] = [];
  for (let i = 0; i < P; i++)
    for (let j = i + 1; j < P; j++) {
      if (focus >= 0 && i !== focus && j !== focus) continue;
      let q = 0;
      for (let t = 0; t < d; t++) q += (X[i][t] - X[j][t]) ** 2;
      if (q < 1) out.push([i, j]);
    }
  return out;
}

// ---------------------------------------------------------------- cubes

export class NCubeCubeModel implements PieceModel {
  readonly target = 1e-20;
  readonly maxStep = 0.25;
  readonly floor: number;
  readonly n: number;
  private readonly maxParts: number;
  private pCell: Int32Array;
  private pCube: Int32Array;
  private loI: Int32Array;
  private loA: Float64Array;
  private loB: Float64Array;
  private hiI: Int32Array;
  private hiA: Float64Array;
  private hiB: Float64Array;
  private loS: Float64Array;
  private hiS: Float64Array;
  private order: Int32Array;

  constructor(
    readonly d: number,
    readonly P: number,
  ) {
    this.n = d - 1;
    // no aligned unit cube fits in a cell of edge < 1; and the volume bound P <= 2d s^n
    this.floor = ncubeBound(d, "cube", P);
    this.maxParts = (this.n + 1) * P;
    const m = this.maxParts * this.n;
    this.pCell = new Int32Array(this.maxParts);
    this.pCube = new Int32Array(this.maxParts);
    this.loI = new Int32Array(m);
    this.loA = new Float64Array(m);
    this.loB = new Float64Array(m);
    this.hiI = new Int32Array(m);
    this.hiA = new Float64Array(m);
    this.hiB = new Float64Array(m);
    this.loS = new Float64Array(m);
    this.hiS = new Float64Array(m);
    this.order = new Int32Array(this.maxParts);
  }

  /**
   * Builds the parts of every cube; returns their count and the line-crossing term (energy, or with `violation` the largest
   * second crossing depth), adding its gradient into g. With the kernel's options the term is weighted and kept for the
   * focus cube only, and its dE/ds (a depth falls by 1/2 per unit of s) is added into ds[0].
   */
  private build(
    st: State,
    s: number,
    g: Float64Array | null,
    violation: boolean,
    opts?: TermOptions,
  ): { count: number; line: number } {
    const { d, n, P } = this;
    const h = s / 2;
    const U = st.U;
    const W = opts?.w ?? null;
    const focus = opts?.focus ?? -1;
    const ds = opts?.ds ?? null;
    let count = 0;
    let line = violation ? -Infinity : 0;
    const dep = new Float64Array(n);
    // the slab's bounds per R^d axis
    const lo4I = new Int32Array(d);
    const lo4A = new Float64Array(d);
    const lo4B = new Float64Array(d);
    const lo4S = new Float64Array(d);
    const hi4I = new Int32Array(d);
    const hi4A = new Float64Array(d);
    const hi4B = new Float64Array(d);
    const hi4S = new Float64Array(d);
    for (let q = 0; q < P; q++) {
      const c = st.cells[q];
      const i = c >> 1;
      const sg = c & 1 ? 1 : -1;
      const ax = axesOf(d, i);
      const o = n * q;
      const wq = W ? W[q] : 1;
      const mine = focus < 0 || q === focus;
      let ncross = 0;
      let kmax = -1;
      for (let k = 0; k < n; k++) {
        dep[k] = Math.abs(U[o + k]) + 0.5 - h;
        if (dep[k] > 0) {
          ncross++;
          if (kmax < 0 || dep[k] > dep[kmax]) kmax = k;
        }
      }
      if (!mine) {
        // another cube's line term is not asked for
      } else if (violation) {
        // the second largest depth: <= 0 exactly when at most one axis crosses
        let a = -Infinity;
        let b = -Infinity;
        for (let k = 0; k < n; k++) {
          if (dep[k] > a) {
            b = a;
            a = dep[k];
          } else if (dep[k] > b) b = dep[k];
        }
        line = Math.max(line, b);
      } else if (ncross >= 2)
        for (let k = 0; k < n; k++)
          if (dep[k] > 0 && k !== kmax) {
            line += wq * dep[k] * dep[k];
            if (g) g[o + k] += 2 * wq * dep[k] * (U[o + k] >= 0 ? 1 : -1);
            if (ds) ds[0] -= wq * dep[k];
          }
      // the part in the anchor cell
      const a = count++;
      this.pCell[a] = c;
      this.pCube[a] = q;
      for (let k = 0; k < n; k++) {
        const x = n * a + k;
        this.loI[x] = o + k;
        this.loA[x] = 1;
        this.loB[x] = -0.5;
        this.loS[x] = 0;
        this.hiI[x] = o + k;
        this.hiA[x] = 1;
        this.hiB[x] = 0.5;
        this.hiS[x] = 0;
        if (dep[k] > 0) {
          if (U[o + k] > 0) {
            this.hiI[x] = -1;
            this.hiA[x] = 0;
            this.hiB[x] = h;
            this.hiS[x] = 0.5;
          } else {
            this.loI[x] = -1;
            this.loA[x] = 0;
            this.loB[x] = -h;
            this.loS[x] = -0.5;
          }
        }
      }
      // a slab beyond each crossed face
      for (let k = 0; k < n; k++) {
        if (!(dep[k] > 0)) continue;
        const tau = U[o + k] > 0 ? 1 : -1;
        const j = ax[k];
        const c2 = 2 * j + (tau > 0 ? 1 : 0);
        const ax2 = axesOf(d, j);
        // along the old normal axis i: from the shared ridge x_i = sg h, depth dd = tau u_k + 1/2 - h
        if (sg > 0) {
          lo4I[i] = o + k;
          lo4A[i] = -tau;
          lo4B[i] = 2 * h - 0.5;
          lo4S[i] = 1;
          hi4I[i] = -1;
          hi4A[i] = 0;
          hi4B[i] = h;
          hi4S[i] = 0.5;
        } else {
          lo4I[i] = -1;
          lo4A[i] = 0;
          lo4B[i] = -h;
          lo4S[i] = -0.5;
          hi4I[i] = o + k;
          hi4A[i] = tau;
          hi4B[i] = 0.5 - 2 * h;
          hi4S[i] = -1;
        }
        for (let m = 0; m < n; m++) {
          if (m === k) continue;
          const gm = ax[m];
          const x = n * a + m;
          lo4I[gm] = this.loI[x];
          lo4A[gm] = this.loA[x];
          lo4B[gm] = this.loB[x];
          lo4S[gm] = this.loS[x];
          hi4I[gm] = this.hiI[x];
          hi4A[gm] = this.hiA[x];
          hi4B[gm] = this.hiB[x];
          hi4S[gm] = this.hiS[x];
        }
        const b = count++;
        this.pCell[b] = c2;
        this.pCube[b] = q;
        for (let p = 0; p < n; p++) {
          const gp = ax2[p];
          const x = n * b + p;
          this.loI[x] = lo4I[gp];
          this.loA[x] = lo4A[gp];
          this.loB[x] = lo4B[gp];
          this.loS[x] = lo4S[gp];
          this.hiI[x] = hi4I[gp];
          this.hiA[x] = hi4A[gp];
          this.hiB[x] = hi4B[gp];
          this.hiS[x] = hi4S[gp];
        }
      }
    }
    return { count, line };
  }

  private bound(I: Int32Array, A: Float64Array, B: Float64Array, x: number, U: Float64Array): number {
    return I[x] < 0 ? B[x] : A[x] * U[I[x]] + B[x];
  }

  private evaluate(st: State, s: number, g: Float64Array | null, violation: boolean, opts?: TermOptions): number {
    const { d, n } = this;
    const ds = opts?.ds ?? null;
    if (ds) ds[0] = 0;
    const { count, line } = this.build(st, s, g, violation, opts);
    const W = opts?.w ?? null;
    const focus = opts?.focus ?? -1;
    const U = st.U;
    let E = line;
    // parts by cell (counting sort)
    const C = 2 * d;
    const first = new Int32Array(C + 1);
    for (let p = 0; p < count; p++) first[this.pCell[p] + 1]++;
    for (let c = 0; c < C; c++) first[c + 1] += first[c];
    const fill = first.slice(0, C);
    for (let p = 0; p < count; p++) this.order[fill[this.pCell[p]]++] = p;
    for (let c = 0; c < C; c++)
      for (let x = first[c]; x < first[c + 1]; x++) {
        const a = this.order[x];
        for (let y = x + 1; y < first[c + 1]; y++) {
          const b = this.order[y];
          if (this.pCube[a] === this.pCube[b]) continue;
          if (focus >= 0 && this.pCube[a] !== focus && this.pCube[b] !== focus) continue;
          let best = Infinity;
          let hiX = -1;
          let loX = -1;
          for (let m = 0; m < n; m++) {
            const ha = this.bound(this.hiI, this.hiA, this.hiB, n * a + m, U);
            const hb = this.bound(this.hiI, this.hiA, this.hiB, n * b + m, U);
            const la = this.bound(this.loI, this.loA, this.loB, n * a + m, U);
            const lb = this.bound(this.loI, this.loA, this.loB, n * b + m, U);
            const ov = Math.min(ha, hb) - Math.max(la, lb);
            if (ov < best) {
              best = ov;
              hiX = ha < hb ? n * a + m : n * b + m;
              loX = la > lb ? n * a + m : n * b + m;
            }
            if (!violation && ov <= 0) break;
          }
          if (violation) {
            E = Math.max(E, best);
            continue;
          }
          if (!(best > 0)) continue;
          const wab = W ? 0.5 * (W[this.pCube[a]] + W[this.pCube[b]]) : 1;
          E += wab * best * best;
          if (ds) ds[0] += 2 * wab * best * (this.hiS[hiX] - this.loS[loX]);
          if (!g) continue;
          if (this.hiI[hiX] >= 0) g[this.hiI[hiX]] += 2 * wab * best * this.hiA[hiX];
          if (this.loI[loX] >= 0) g[this.loI[loX]] -= 2 * wab * best * this.loA[loX];
        }
      }
    return E;
  }

  energy(st: State, s: number, g: Float64Array): number {
    g.fill(0);
    return this.evaluate(st, s, g, false);
  }

  terms(st: State, s: number, g: Float64Array | null, violation: boolean, o: TermOptions): number {
    if (g) g.fill(0);
    return this.evaluate(st, s, g, violation, o);
  }

  violation(st: State, s: number): number {
    return this.evaluate(st, s, null, true);
  }

  startEdge(): number {
    // (2d k^n aligned cubes fit at s = k)
    return Math.max(2, Math.ceil(Math.pow(this.P / (2 * this.d), 1 / this.n)));
  }

  randomState(s: number, rand: () => number): State {
    const { d, n, P } = this;
    const cells = new Int32Array(P);
    const U = new Float64Array(n * P);
    // centres within 1/2 of the cell's centre
    const shrink = Math.min(1, 1 / s);
    for (let q = 0; q < P; q++) {
      cells[q] = randInt(rand, 2 * d);
      U.set(randomPoint(n, s, rand, shrink), n * q);
    }
    return { cells, U };
  }

  hop(st: State, s: number, rand: () => number): State {
    const { d, n, P } = this;
    const out: State = { cells: Int32Array.from(st.cells), U: Float64Array.from(st.U) };
    const k = 1 + randInt(rand, 3);
    for (let t = 0; t < k; t++) {
      const q = randInt(rand, P);
      out.cells[q] = randInt(rand, 2 * d);
      out.U.set(randomPoint(n, s, rand, 1), n * q);
    }
    return out;
  }

  constructedStart(): null {
    return null;
  }

  stateOf(config: BoundaryConfig): State | null {
    return plainState(this.d, config, this.P);
  }

  pieces(st: State, s: number): BoundaryPiece[] {
    return clampedPieces(this.d, st, s);
  }
}
