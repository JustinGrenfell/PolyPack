/**
 * The nD kernel for the generic search methods (packages/methods/src/kernel.ts): P balls (diameter 1) on the boundary of the
 * simplex, orthoplex or demicube, d = 5..10, the twin of ncube/kernel.ts with the chain engine (chains.ts) in place of the lattice words.
 *
 * State. aux holds each piece's anchor: a cell and a centre in that cell's local frame, recorded at the edge s0 (aux.s0). At the state's
 * size the anchor is that centre scaled by size / s0. x holds each piece's displacement from its anchor, n = d - 1 lengths in the anchor's
 * frame; the piece's centre is the end of the straight segment from the anchor along it, unfolded across every ridge it crosses
 * (move.ts walkND). rebase moves every anchor to its piece's centre and zeroes x. The size is the edge s.
 *
 * Energy: the model's penalty energy at the centres (model.ts), weighted per piece. Its gradient in a piece's x is M^T times its
 * gradient at the centre, M the linear part of the reflections crossed. In the size it is the model's dE/ds at fixed local
 * coordinates plus the centres' own motion: the geometry at edge s is the edge-1 geometry scaled by s, so a centre y = s G(a, v/s) moves
 * as dy/ds = (y - M v)/s (ambient), and the local coordinates u = f . (y - s o) as f . (dy/ds - o). A segment that meets a hinge
 * has no continuation: such a displacement (a set of measure zero) leaves its piece at the anchor. So does one longer than
 * MAX_WALK_EDGES edges.
 *
 * finish: the search's finalize stage (finalizeState) with `shrink`; validate: the server's rule (rule.ts).
 * Capabilities: "seeds" (the exact layouts and the sequential layout of seeds.ts) and "shape", a no-op (a ball is its own inscribed ball).
 */
import type { Finished, KernelCapability, KernelState, Rng, SpaceKernel, Validation } from "@polypack/methods";
import { finalizeState } from "../search/search";
import type { State } from "../search/model";
import type { BoundaryConfig } from "../types";
import { BALL_RADIUS } from "../units";
import { getND } from "./container";
import { NDBallModel } from "./model";
import { walkND } from "./move";
import { ND_MAX_P, type NDName, isNDName } from "./names";
import { ndFormErrors, ndServerRule } from "./rule";
import { ndExactLayout, ndSequentialLayout } from "./seeds";
import { type NDProblem, validateND } from "./validate";

/** Energy evaluations one finish() and one validate() cost (measured, BOUNDARY.md 4D.8). */
export const ND_FINISH_COST = 500;
export const ND_VALIDATE_COST = 300;

interface Anchors {
  cells: Int32Array;
  U: Float64Array;
  s0: number;
}

/** A displacement longer than this many edges leaves its piece at the anchor. */
const MAX_WALK_EDGES = 64;
const POSE_TRIES = 16;

export function ndKernel(problem: { polytope: NDName; piece: "ball"; P: number }): SpaceKernel {
  if (
    !problem ||
    !isNDName(problem.polytope) ||
    problem.piece !== "ball" ||
    !Number.isInteger(problem.P) ||
    problem.P < 1
  )
    throw new Error("ndKernel: invalid problem");
  if (problem.P > ND_MAX_P)
    throw new Error(`ndKernel: P = ${problem.P} is above the largest of the ${problem.polytope}, ${ND_MAX_P}`);
  const C = getND(problem.polytope);
  const n = C.n;
  const D = C.D;
  const prob: NDProblem = { polytope: problem.polytope, piece: "ball", P: problem.P };
  const { P } = prob;
  const model = new NDBallModel(prob.polytope, P);
  const nv = n * P;
  const step = BALL_RADIUS;
  const placed: State = { cells: new Int32Array(P), U: new Float64Array(nv) };
  const Ms = new Float64Array(P * D * D);
  const ends = Array.from({ length: P }, () => new Float64Array(D));
  const mv = Array.from({ length: P }, () => new Float64Array(D));
  const gU = new Float64Array(nv);
  const ds = new Float64Array(1);
  const stuck = new Uint8Array(P);
  const bound = model.floor;
  const anchorsOf = (s: KernelState) => s.aux as Anchors;

  /** The ambient vector of the displacement x[o..o+n) in the frame of cell c. */
  const ambientVec = (c: number, x: ArrayLike<number>, o: number, out: Float64Array): Float64Array => {
    const fr = C.cell(c).frame;
    out.fill(0);
    for (let k = 0; k < n; k++) for (let r = 0; r < D; r++) out[r] += x[o + k] * fr[k][r];
    return out;
  };

  const place = (s: KernelState): State => {
    const a = anchorsOf(s);
    const size = s.size;
    const k = size === a.s0 ? 1 : size / a.s0;
    const v = new Float64Array(D);
    const start = new Float64Array(n);
    for (let q = 0; q < P; q++) {
      const o = n * q;
      for (let t = 0; t < n; t++) start[t] = a.U[o + t] * k;
      const xs = C.toAmbient(a.cells[q], start, size);
      ambientVec(a.cells[q], s.x, o, v);
      let len2 = 0;
      for (let r = 0; r < D; r++) len2 += v[r] * v[r];
      const M = Ms.subarray(q * D * D, (q + 1) * D * D);
      let w = null as ReturnType<typeof walkND> | null;
      if (len2 <= (MAX_WALK_EDGES * size) ** 2) w = walkND(C, size, a.cells[q], xs, v, M);
      if (w && !w.hit) {
        stuck[q] = 0;
        placed.cells[q] = w.cell;
        C.toLocal(w.cell, w.x, size, placed.U.subarray(o, o + n));
        ends[q].set(w.x);
        mv[q].fill(0);
        for (let r = 0; r < D; r++) {
          let t = 0;
          for (let c = 0; c < D; c++) t += M[r * D + c] * v[c];
          mv[q][r] = t;
        }
      } else {
        stuck[q] = 1;
        placed.cells[q] = a.cells[q];
        for (let t = 0; t < n; t++) placed.U[o + t] = start[t];
        ends[q].set(xs);
        mv[q].fill(0);
        M.fill(0);
        for (let r = 0; r < D; r++) M[r * D + r] = 1;
      }
    }
    return placed;
  };

  const stateOf = (st: State, size: number): KernelState => ({
    x: new Float64Array(nv),
    size,
    aux: { cells: Int32Array.from(st.cells), U: Float64Array.from(st.U), s0: size } satisfies Anchors,
  });

  const accept = (config: BoundaryConfig) => validateND(prob, config, { tol: 0 }).ok;
  const capabilities: KernelCapability[] = ["shape", "seeds"];
  const poseV = new Float64Array(D);

  const kernel: SpaceKernel = {
    space: "boundary",
    id: `boundary/${prob.polytope}/ball/${P}`,
    P,
    dim: n,
    steps: new Float64Array(n).fill(step),
    inradius: step,
    globals: 0,
    globalSteps: new Float64Array(0),
    bound,
    finishCost: ND_FINISH_COST,
    validateCost: ND_VALIDATE_COST,
    capabilities,

    fromConfig(config) {
      if (ndFormErrors(prob, config).length) return null;
      const st = model.stateOf(config as BoundaryConfig);
      return st ? stateOf(st, (config as BoundaryConfig).s) : null;
    },

    toConfig(s): BoundaryConfig {
      return { s: s.size, pieces: model.pieces(place(s), s.size) };
    },

    clone(s) {
      const a = anchorsOf(s);
      return {
        x: Float64Array.from(s.x),
        size: s.size,
        aux: { cells: Int32Array.from(a.cells), U: Float64Array.from(a.U), s0: a.s0 } satisfies Anchors,
      };
    },

    rebase(s) {
      const st = place(s);
      const a = anchorsOf(s);
      a.cells.set(st.cells);
      a.U.set(st.U);
      a.s0 = s.size;
      s.x.fill(0);
    },

    randomStart(rng) {
      const s = model.startEdge();
      return stateOf(model.randomState(s, rng), s);
    },

    randomPose(s, i, rng) {
      // the end of a random geodesic: a direction uniform on the sphere, a length up to 1.5 edges with the density of a ball, tried by
      // walking it; a walk that meets a hinge is drawn again, and after POSE_TRIES misses the piece stays where it is
      const a = anchorsOf(s);
      const size = s.size;
      const k = size === a.s0 ? 1 : size / a.s0;
      const from = new Float64Array(n);
      for (let t = 0; t < n; t++) from[t] = a.U[n * i + t] * k;
      const xs = C.toAmbient(a.cells[i], from, size);
      const v = new Float64Array(n);
      for (let t = 0; t < POSE_TRIES; t++) {
        let r2 = 0;
        for (let q = 0; q < n; q++) {
          const g = Math.sqrt(-2 * Math.log(1 - rng())) * Math.cos(2 * Math.PI * rng());
          v[q] = g;
          r2 += g * g;
        }
        const L = (1.5 * size * Math.pow(rng(), 1 / n)) / Math.sqrt(Math.max(r2, 1e-12));
        for (let q = 0; q < n; q++) v[q] *= L;
        ambientVec(a.cells[i], v, 0, poseV);
        if (walkND(C, size, a.cells[i], xs, poseV).hit) continue;
        for (let q = 0; q < n; q++) s.x[n * i + q] = v[q];
        return;
      }
    },

    jiggle(s, i, amount, rng) {
      for (let m = 0; m < n; m++) s.x[n * i + m] += (2 * rng() - 1) * amount * step;
    },

    setSize(s, size) {
      const k = size / s.size;
      for (let i = 0; i < nv; i++) s.x[i] *= k;
      s.size = size;
    },

    energy(s, grad, weights) {
      const st = place(s);
      const size = s.size;
      if (!grad) return model.terms(st, size, null, false, { w: weights ?? null });
      const e = model.terms(st, size, gU, false, { w: weights ?? null, ds });
      let dsize = ds[0];
      const G = new Float64Array(D);
      const G0 = new Float64Array(D);
      for (let q = 0; q < P; q++) {
        const o = n * q;
        const cEnd = st.cells[q];
        const frEnd = C.cell(cEnd).frame;
        G.fill(0);
        for (let k = 0; k < n; k++) for (let r = 0; r < D; r++) G[r] += gU[o + k] * frEnd[k][r];
        // d y / d s (ambient) = (y - M v) / size, and the local coordinates lose the centre's s o
        const oc = C.cell(cEnd).o;
        let t = 0;
        for (let r = 0; r < D; r++) t += G[r] * ((ends[q][r] - mv[q][r]) / size - oc[r]);
        dsize += t;
        const M = Ms.subarray(q * D * D, (q + 1) * D * D);
        // dE/dv = M^T G, then onto the anchor's frame
        for (let r = 0; r < D; r++) {
          let u = 0;
          for (let c = 0; c < D; c++) u += M[c * D + r] * G[c];
          G0[r] = u;
        }
        const frA = C.cell(anchorsOf(s).cells[q]).frame;
        for (let k = 0; k < n; k++) {
          let u = 0;
          for (let r = 0; r < D; r++) u += frA[k][r] * G0[r];
          grad[o + k] = stuck[q] ? 0 : u;
        }
      }
      grad[nv] = dsize;
      return e;
    },

    pieceEnergy(s, i, weights) {
      return model.terms(place(s), s.size, null, false, { w: weights ?? null, focus: i });
    },

    clear(s, i, tol) {
      return model.terms(place(s), s.size, null, true, { focus: i }) <= tol;
    },

    finish(s): Finished | null {
      const st = place(s);
      const copy: State = { cells: Int32Array.from(st.cells), U: Float64Array.from(st.U) };
      const f = finalizeState(model, copy, s.size, accept, { shrink: true });
      return f ? { config: f.config, score: f.score } : null;
    },

    validate(config): Validation {
      return ndServerRule(prob, config, bound);
    },

    seeds(rng: Rng, count: number) {
      const out: KernelState[] = [];
      for (let t = 0; t < count; t++) {
        // the exact layout (when one is known) and the sequential layout alternate, the exact one first
        const built =
          (t % 2 === 0 ? ndExactLayout(prob.polytope, P, rng) : null) ?? ndSequentialLayout(prob.polytope, P, rng);
        out.push(stateOf(built.st, built.s));
      }
      return out;
    },
  };
  kernel.setShape = () => {};
  return kernel;
}
