/**
 * The d-cube's kernel for the generic search methods (packages/methods/src/kernel.ts), d = 5..10: P balls (diameter 1) or
 * axis-aligned unit cubes on the boundary of the d-cube, the d-cube twin of ../kernel.ts with n = d - 1 variables per piece.
 *
 * State. aux holds each piece's anchor: a cell and a centre in that cell's local frame, recorded at the edge s0 (aux.s0). At
 * the state's size the anchor is that centre scaled by size / s0. x holds each piece's displacement from its anchor, n lengths
 * in the anchor's frame; the piece's centre is the end of the straight segment from the anchor along it, unfolded across every
 * face it crosses (move.ts walkN). rebase moves every anchor to its piece's centre and zeroes x. The size is the edge s.
 *
 * Energy: the model's penalty energy at the centres (model.ts), weighted per piece. Its gradient in a piece's x is M^T times
 * its gradient at the centre, M the signed permutation of the gluings crossed. In the size it is the model's dE/ds at fixed
 * centres plus the centres' own motion: the geometry at edge s is the edge-1 geometry scaled by s, so a centre y = s G(a, v/s)
 * moves as dy/ds = (y - M v) / s. A segment that meets a hinge has no continuation: such a displacement (a set of measure
 * zero) leaves its piece at the anchor.
 *
 * finish: the search's finalize stage (finalizeState) with `shrink`; validate: the server's rule (rule.ts).
 *
 * Capabilities: "seeds" (the cell layouts of seeds.ts; for cubes also the ring cut; for balls the ring cut when the validator
 * accepts it), and for balls "shape", a no-op.
 */
import type { Finished, KernelCapability, KernelState, Rng, SpaceKernel, Validation } from "@polypack/methods";
import { finalizeState } from "../search/search";
import { randInt } from "../search/random";
import type { PieceModel, State } from "../search/model";
import type { BoundaryConfig } from "../types";
import { BALL_RADIUS } from "../units";
import { NCUBE_MAX_P } from "./bounds";
import { NCubeName, isNCubeName, ncubeDim } from "./lattice";
import { NCubeBallModel, NCubeCubeModel, randomPoint } from "./model";
import { walkN } from "./move";
import { ncubeFormErrors, ncubeServerRule } from "./rule";
import { cellLayout, ringLayout } from "./seeds";
import { validateNCube } from "./validate";

/** Energy evaluations one finish() and one validate() cost (a validator call plus 20 to 45 evaluations of the model's rule). */
export const NCUBE_FINISH_COST = 1000;
export const NCUBE_VALIDATE_COST = 1000;

interface Anchors {
  cells: Int32Array;
  U: Float64Array;
  s0: number;
}

/** A displacement longer than this many edges leaves its piece at the anchor. */
const MAX_WALK_EDGES = 64;
const POSE_TRIES = 16;

export function ncubeKernel(problem: { polytope: NCubeName; piece: "ball" | "cube"; P: number }): SpaceKernel {
  if (
    !problem ||
    !isNCubeName(problem.polytope) ||
    (problem.piece !== "ball" && problem.piece !== "cube") ||
    !Number.isInteger(problem.P) ||
    problem.P < 1
  )
    throw new Error("ncubeKernel: invalid problem");
  if (problem.P > NCUBE_MAX_P[problem.polytope])
    throw new Error(
      `ncubeKernel: P = ${problem.P} is above the largest of the ${problem.polytope}, ${NCUBE_MAX_P[problem.polytope]}`,
    );
  const d = ncubeDim(problem.polytope);
  const n = d - 1;
  const C = 2 * d;
  const prob = { polytope: problem.polytope, piece: problem.piece, P: problem.P };
  const { P, piece } = prob;
  const ball = piece === "ball";
  const model: PieceModel = ball ? new NCubeBallModel(d, P) : new NCubeCubeModel(d, P);
  const nv = n * P;
  const step = ball ? BALL_RADIUS : 0.5;
  const placed: State = { cells: new Int32Array(P), U: new Float64Array(nv) };
  const perm = new Int32Array(nv);
  const sgn = new Float64Array(nv);
  const gU = new Float64Array(nv);
  const ds = new Float64Array(1);
  const stuck = new Uint8Array(P);
  const poseU = new Float64Array(n);
  const start = new Float64Array(n);
  const bound = model.floor;

  const anchorsOf = (s: KernelState) => s.aux as Anchors;

  /** The pieces of s at their centres (scratch `placed`), with each map's signed permutation in perm/sgn. */
  const place = (s: KernelState): State => {
    const a = anchorsOf(s);
    const size = s.size;
    const k = size === a.s0 ? 1 : size / a.s0;
    for (let q = 0; q < P; q++) {
      const o = n * q;
      for (let t = 0; t < n; t++) start[t] = a.U[o + t] * k;
      const c = walkLong(d, a.cells[q], start, s.x, o, size, placed.U, o, perm, sgn, o);
      stuck[q] = c >= 0 ? 0 : 1;
      if (c >= 0) placed.cells[q] = c;
      else {
        placed.cells[q] = a.cells[q];
        for (let t = 0; t < n; t++) {
          placed.U[o + t] = start[t];
          perm[o + t] = t;
          sgn[o + t] = 0;
        }
      }
    }
    return placed;
  };

  const stateOf = (st: State, size: number): KernelState => ({
    x: new Float64Array(nv),
    size,
    aux: { cells: Int32Array.from(st.cells), U: Float64Array.from(st.U), s0: size } satisfies Anchors,
  });

  const accept = (config: BoundaryConfig) => validateNCube(prob, config, { tol: 0 }).ok;
  const capabilities: KernelCapability[] = ball ? ["shape", "seeds"] : ["seeds"];

  const kernel: SpaceKernel = {
    space: "boundary",
    id: `boundary/${prob.polytope}/${piece}/${P}`,
    P,
    dim: n,
    steps: new Float64Array(n).fill(step),
    inradius: step,
    globals: 0,
    globalSteps: new Float64Array(0),
    bound,
    finishCost: NCUBE_FINISH_COST,
    validateCost: NCUBE_VALIDATE_COST,
    capabilities,

    fromConfig(config) {
      if (ncubeFormErrors(prob, config).length) return null;
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
      if (ball) {
        const s = model.startEdge();
        return stateOf(model.randomState(s, rng), s);
      }
      // cubes: an edge at which the cubes fill at most a quarter of the boundary, each cube still inside its cell
      const s = Math.max(model.startEdge(), Math.pow((2 * P) / d, 1 / n));
      const cells = new Int32Array(P);
      const U = new Float64Array(nv);
      for (let q = 0; q < P; q++) {
        cells[q] = randInt(rng, C);
        U.set(randomPoint(n, s, rng, 1 - 1 / s), n * q);
      }
      return stateOf({ cells, U }, s);
    },

    randomPose(s, i, rng) {
      // the end of a random geodesic: a direction uniform on the sphere, a length up to 2.5 edges with the density of a ball,
      // tried by walking it; a walk that meets a hinge is drawn again, and after POSE_TRIES misses the piece stays
      const a = anchorsOf(s);
      const size = s.size;
      const k = size === a.s0 ? 1 : size / a.s0;
      const from = new Float64Array(n);
      for (let t = 0; t < n; t++) from[t] = a.U[n * i + t] * k;
      const v = new Float64Array(n);
      for (let t = 0; t < POSE_TRIES; t++) {
        let r2 = 0;
        for (let q = 0; q < n; q++) {
          // a Gaussian direction (Box-Muller from the stream)
          const g = Math.sqrt(-2 * Math.log(1 - rng())) * Math.cos(2 * Math.PI * rng());
          v[q] = g;
          r2 += g * g;
        }
        const L = (2.5 * size * Math.pow(rng(), 1 / n)) / Math.sqrt(Math.max(r2, 1e-12));
        for (let q = 0; q < n; q++) v[q] *= L;
        if (walkLong(d, a.cells[i], from, v, 0, size, poseU, 0, null, null, 0) < 0) continue;
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
      const x = s.x;
      for (let q = 0; q < P; q++) {
        const o = n * q;
        const free = !stuck[q];
        for (let r = 0; r < n; r++) {
          // dE/dv = M^T dE/dy: v_end[r] = sgn[r] v_start[perm[r]]
          if (free) grad[o + perm[o + r]] = 0;
          else grad[o + r] = 0;
        }
        for (let r = 0; r < n; r++) {
          if (free) grad[o + perm[o + r]] += sgn[o + r] * gU[o + r];
          // the centre's motion with the size: dy/ds = (y - M v) / s (y / s for a piece at its anchor)
          const mv = free ? sgn[o + r] * x[o + perm[o + r]] : 0;
          dsize += (gU[o + r] * (st.U[o + r] - mv)) / size;
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
      return ncubeServerRule(prob, config, bound);
    },

    seeds(rng: Rng, count: number) {
      const out: KernelState[] = [];
      for (let t = 0; t < count; t++) {
        // the ring cut and the cell layout alternate (the ring first); a ball ring is kept only when the validator accepts it
        let built: { st: State; s: number } | null = null;
        if (t % 2 === 0) {
          built = ringLayout(d, piece, P, rng);
          if (built && ball && !accept({ s: built.s, pieces: model.pieces(built.st, built.s) })) built = null;
        }
        built ??= cellLayout(d, piece, P, rng);
        out.push(stateOf(built.st, built.s));
      }
      return out;
    },
  };
  if (ball) kernel.setShape = () => {};
  return kernel;
}

/**
 * walkN for a displacement of any length: one call when it is shorter than 8 edges, else in equal parts, each started
 * where and in the frame the previous one ended (the parts compose; their signed permutations multiply). Returns the end
 * cell, or -1 when the segment meets a hinge.
 */
function walkLong(
  d: number,
  cell: number,
  p: Float64Array,
  v: ArrayLike<number>,
  vo: number,
  s: number,
  outU: Float64Array,
  ou: number,
  perm: Int32Array | null,
  sgn: Float64Array | null,
  po: number,
): number {
  const n = d - 1;
  let len2 = 0;
  for (let k = 0; k < n; k++) len2 += v[vo + k] * v[vo + k];
  const len = Math.sqrt(len2);
  if (!(len <= MAX_WALK_EDGES * s)) return -1;
  const parts = Math.ceil(len / (8 * s));
  if (!(parts > 1)) return walkN(d, cell, p, 0, v, vo, s, outU, ou, perm, sgn, po);
  let dvec = Float64Array.from({ length: n }, (_, k) => v[vo + k] / parts);
  const cur = Float64Array.from(p);
  // total signed permutation: v_end[r] = tsg[r] v_start[tp[r]]
  const tp = Int32Array.from({ length: n }, (_, k) => k);
  const ts = new Float64Array(n).fill(1);
  const pp = new Int32Array(n);
  const ps = new Float64Array(n);
  let c = cell;
  for (let t = 0; t < parts; t++) {
    c = walkN(d, c, cur, 0, dvec, 0, s, cur, 0, pp, ps, 0);
    if (c < 0) return -1;
    const next = new Float64Array(n);
    const np = new Int32Array(n);
    const nsg = new Float64Array(n);
    for (let r = 0; r < n; r++) {
      next[r] = ps[r] * dvec[pp[r]];
      np[r] = tp[pp[r]];
      nsg[r] = ps[r] * ts[pp[r]];
    }
    dvec = next;
    tp.set(np);
    ts.set(nsg);
  }
  for (let k = 0; k < n; k++) outU[ou + k] = cur[k];
  if (perm && sgn)
    for (let k = 0; k < n; k++) {
      perm[po + k] = tp[k];
      sgn[po + k] = ts[k];
    }
  return c;
}
