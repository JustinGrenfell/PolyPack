/**
 * The boundary kernel (kernel.ts) for the Platonic solids: P rotated solids on the boundary of a regular 4-polytope, the tesseract's
 * axis-aligned cube excepted (that one is kernel.ts's). docs/spaces/BOUNDARY.md section 4C and 14.
 *
 * State. As kernel.ts, with an orientation: aux holds each piece's anchor (a cell, a centre recorded at the edge s0, and an orientation R_a,
 * an orthogonal 3 x 3 matrix in the anchor's frame). x has dim = 6 numbers per piece: the displacement v (3 lengths, as kernel.ts) and a
 * rotation vector w (3 radians) about the anchor's orientation. The piece's centre is where the straight path from the anchor along v ends,
 * unfolded across every ridge it crosses (the chart's exponential map); its orientation there is  M exp(w) R_a,  M the linear part of the
 * gluings crossed (a piece left at its anchor because its path met a hinge keeps the rotation, with M = I). rebase moves the anchors to the
 * pieces and zeroes x. The energy is the model's (search/solids.ts); its gradient in w is  Jl(w)^T det(M) M^T  times the model's gradient in
 * the rotation increment at the end point (Jl the left Jacobian of the exponential; det M = -1 when a gluing reverses handedness).
 *
 * Moves change a piece's x only, so writing x back undoes them. randomPose draws a target cell and point (as kernel.ts) and a uniform
 * orientation; w is then the rotation that gives that orientation after the walk. swap gives piece i the orientation of piece j and vice
 * versa: every Platonic solid is mirror-symmetric, so an orientation that would need a handedness change is replaced by its mirror image
 * (the same piece), which always lies in the same coset as the walk requires.
 *
 * Capabilities: "seeds" (the cells themselves when the solid is the container's cell and P is at most the cell count; else a sequential
 * random layout accepted piece by piece by the model's violation) and "swap". No "shape": the model has no soft solids.
 */
import type { Finished, KernelCapability, KernelState, Rng, SpaceKernel, Validation } from "@polypack/methods";
import { pieceBounds } from "./bounds";
import { formatBoundaryId } from "./ids";
import { displacementTo } from "./kernel-pose";
import { walk } from "./kernel";
import { BOUNDARY_SERVER_RULE, serverForm, serverRule } from "./kernel-validate";
import { getPolytope } from "./polytope";
import { type Charts, chartsOf, randomCellPoint } from "./search/chart";
import type { State } from "./search/model";
import { randInt } from "./search/random";
import { expRot, leftJacobian, logRot, mul9, det9, orthonormalize, randomRotation } from "./search/rot";
import { finalizeState } from "./search/search";
import { SolidModel } from "./search/solids";
import type { BoundaryConfig, BoundaryProblem, PolytopeName, PolytopePieceName } from "./types";
import { validateBoundary } from "./validate";

interface Anchors {
  cells: Int32Array;
  U: Float64Array;
  /** 9 numbers per piece: the anchor orientation. */
  R: Float64Array;
  s0: number;
}

const POSE_TRIES = 16;

/**
 * Energy evaluations one finish() and one validate() cost, in energy() calls with their gradients. [CALC] Process CPU time on an Apple Silicon
 * laptop, on squeezed random starts of 4 and 12 solids of the five kinds on the six containers: a finish (a validator call and 20 to 45 evaluations
 * of the model's rule) took 12 to 55 energy() calls, a validator call 0.1 to 1 on random states and up to about 15 on a crowd of touching solids.
 */
export const SOLID_FINISH_COST = 60;
export const SOLID_VALIDATE_COST = 15;

export function solidKernel(prob: BoundaryProblem & { piece: PolytopePieceName }): SpaceKernel {
  const { P, piece } = prob;
  const poly = getPolytope(prob.polytope as PolytopeName);
  const charts: Charts = chartsOf(poly);
  const model = new SolidModel(poly, piece, P);
  const shape = model.shape;
  const dim = 6;
  const n = dim * P;
  // a step: the solid's inradius for a displacement, and the turn that moves its vertices about as far (inradius / circumradius)
  const step = shape.inradius;
  const rotStep = shape.inradius / shape.circumradius;

  const placed: State = { cells: new Int32Array(P), U: new Float64Array(3 * P), R: new Float64Array(9 * P) };
  const Ms = new Float64Array(9 * P);
  const gU = new Float64Array(6 * P);
  const ds = new Float64Array(1);
  const start = new Float64Array(3);
  const partM = new Float64Array(9);
  const stuck = new Uint8Array(P);
  const poseU = new Float64Array(3);
  const poseM = new Float64Array(9);
  const eRot = new Float64Array(9);
  const tmp = new Float64Array(9);
  const jac = new Float64Array(9);
  const serverBound = pieceBounds(prob.polytope as PolytopeName, piece, P).value;
  const anchorsOf = (s: KernelState) => s.aux as Anchors;

  /** The pieces of s at their centres, with their orientations; each map's linear part in Ms (identity for a piece left at its anchor). */
  const place = (s: KernelState): State => {
    const a = anchorsOf(s);
    const size = s.size;
    const k = size === a.s0 ? 1 : size / a.s0;
    const x = s.x;
    for (let q = 0; q < P; q++) {
      const o = 3 * q;
      const xo = dim * q;
      start[0] = a.U[o] * k;
      start[1] = a.U[o + 1] * k;
      start[2] = a.U[o + 2] * k;
      const ax = start[0], ay = start[1], az = start[2];
      const c = walk(charts, a.cells[q], start, x[xo], x[xo + 1], x[xo + 2], size, placed.U, o, Ms, 9 * q, partM);
      stuck[q] = c >= 0 ? 0 : 1;
      if (c >= 0) placed.cells[q] = c;
      else {
        placed.cells[q] = a.cells[q];
        placed.U[o] = ax;
        placed.U[o + 1] = ay;
        placed.U[o + 2] = az;
        Ms.set([1, 0, 0, 0, 1, 0, 0, 0, 1], 9 * q);
      }
      expRot(x[xo + 3], x[xo + 4], x[xo + 5], eRot);
      mul9(eRot, 0, a.R, 9 * q, tmp, 0);
      mul9(Ms, 9 * q, tmp, 0, placed.R!, 9 * q);
      orthonormalize(placed.R!, 9 * q);
    }
    return placed;
  };

  const stateOf = (st: State, size: number): KernelState => ({
    x: new Float64Array(n),
    size,
    aux: { cells: Int32Array.from(st.cells), U: Float64Array.from(st.U), R: Float64Array.from(st.R!), s0: size } satisfies Anchors,
  });

  const accept = (config: BoundaryConfig) => validateBoundary(prob, config, { tol: 0 }).ok;
  const capabilities: KernelCapability[] = ["seeds", "swap"];

  /** A uniformly random orientation reached from the anchor orientation Ra by the map M: w with M exp(w) Ra = a rotation of the solid. */
  const rotationVectorFor = (M: Float64Array, mo: number, Ra: Float64Array, ro: number, target: number[]): [number, number, number] => {
    // exp(w) = M^T T Ra^T; it must be proper, so a target of the wrong handedness is replaced by its mirror image (the same piece)
    const T = target.slice();
    const mT = [M[mo], M[mo + 3], M[mo + 6], M[mo + 1], M[mo + 4], M[mo + 7], M[mo + 2], M[mo + 5], M[mo + 8]];
    const RaT = [Ra[ro], Ra[ro + 3], Ra[ro + 6], Ra[ro + 1], Ra[ro + 4], Ra[ro + 7], Ra[ro + 2], Ra[ro + 5], Ra[ro + 8]];
    const t1 = new Array<number>(9);
    const t2 = new Array<number>(9);
    mul9(mT, 0, T, 0, t1, 0);
    mul9(t1, 0, RaT, 0, t2, 0);
    if (det9(t2) < 0) {
      mul9(T, 0, shape.mirror, 0, t1, 0);
      mul9(mT, 0, t1, 0, t2, 0);
      mul9(t2, 0, RaT, 0, t1, 0);
      for (let k = 0; k < 9; k++) t2[k] = t1[k];
    }
    orthonormalize(t2);
    return logRot(t2);
  };

  const kernel: SpaceKernel = {
    space: "boundary",
    id: formatBoundaryId(prob),
    P,
    dim,
    steps: Float64Array.from([step, step, step, rotStep, rotStep, rotStep]),
    inradius: step,
    globals: 0,
    globalSteps: new Float64Array(0),
    // the proven bound max(volume, one piece) of docs/spaces/BOUNDARY.md 4C.4
    bound: model.floor,
    finishCost: SOLID_FINISH_COST,
    validateCost: SOLID_VALIDATE_COST,
    capabilities,

    fromConfig(config) {
      if (serverForm(prob, config).length) return null;
      const c = config as BoundaryConfig;
      const st = model.stateOf(c);
      return st ? stateOf(st, c.s) : null;
    },

    toConfig(s): BoundaryConfig {
      return { s: s.size, pieces: model.pieces(place(s), s.size) };
    },

    clone(s) {
      const a = anchorsOf(s);
      return {
        x: Float64Array.from(s.x),
        size: s.size,
        aux: { cells: Int32Array.from(a.cells), U: Float64Array.from(a.U), R: Float64Array.from(a.R), s0: a.s0 } satisfies Anchors,
      };
    },

    rebase(s) {
      const st = place(s);
      const a = anchorsOf(s);
      a.cells.set(st.cells);
      a.U.set(st.U);
      a.R.set(st.R!);
      a.s0 = s.size;
      s.x.fill(0);
    },

    randomStart(rng) {
      const s = model.startEdge();
      return stateOf(model.randomState(s, rng), s);
    },

    randomPose(s, i, rng) {
      const a = anchorsOf(s);
      const size = s.size;
      const k = size === a.s0 ? 1 : size / a.s0;
      const from = [a.U[3 * i] * k, a.U[3 * i + 1] * k, a.U[3 * i + 2] * k];
      const tol = 1e-9 * Math.max(1, size);
      const target = randomRotation(rng);
      for (let t = 0; t < POSE_TRIES; t++) {
        const cell = randInt(rng, poly.cells.length);
        const u = randomCellPoint(poly, cell, size, rng, 1);
        const v = displacementTo(a.cells[i], from, cell, u, size, poly);
        if (!v) continue;
        start.set(from);
        const c = walk(charts, a.cells[i], start, v[0], v[1], v[2], size, poseU, 0, poseM, 0, partM);
        if (c !== cell || !(Math.max(Math.abs(poseU[0] - u[0]), Math.abs(poseU[1] - u[1]), Math.abs(poseU[2] - u[2])) <= tol)) continue;
        s.x.set(v, dim * i);
        s.x.set(rotationVectorFor(poseM, 0, a.R, 9 * i, target), dim * i + 3);
        return;
      }
      if (poly.name === "tesseract") return;
      for (let t = 0; t < POSE_TRIES; t++) {
        let d0 = 0, d1 = 0, d2 = 0, r2 = 0;
        do {
          d0 = 2 * rng() - 1;
          d1 = 2 * rng() - 1;
          d2 = 2 * rng() - 1;
          r2 = d0 * d0 + d1 * d1 + d2 * d2;
        } while (r2 > 1 || r2 < 1e-6);
        const L = (2.5 * size * Math.cbrt(rng())) / Math.sqrt(r2);
        start.set(from);
        if (walk(charts, a.cells[i], start, d0 * L, d1 * L, d2 * L, size, poseU, 0, poseM, 0, partM) < 0) continue;
        s.x[dim * i] = d0 * L;
        s.x[dim * i + 1] = d1 * L;
        s.x[dim * i + 2] = d2 * L;
        s.x.set(rotationVectorFor(poseM, 0, a.R, 9 * i, target), dim * i + 3);
        return;
      }
    },

    jiggle(s, i, amount, rng) {
      for (let m = 0; m < 3; m++) s.x[dim * i + m] += (2 * rng() - 1) * amount * step;
      for (let m = 3; m < 6; m++) s.x[dim * i + m] += (2 * rng() - 1) * amount * rotStep;
    },

    setSize(s, size) {
      // the anchors scale with the size; scaling the displacements too scales every centre exactly; rotations are scale-free
      const k = size / s.size;
      for (let q = 0; q < P; q++) for (let m = 0; m < 3; m++) s.x[dim * q + m] *= k;
      s.size = size;
    },

    energy(s, grad, weights) {
      const st = place(s);
      const size = s.size;
      if (!grad) return model.terms(st, size, null, false, { w: weights ?? null });
      const e = model.terms(st, size, gU, false, { w: weights ?? null, ds });
      let dsize = ds[0];
      const x = s.x;
      const off = 3 * P;
      for (let q = 0; q < P; q++) {
        const o = 3 * q;
        const xo = dim * q;
        const m = 9 * q;
        const ga = gU[o], gb = gU[o + 1], gc = gU[o + 2];
        // dE/dv = M^T dE/dy
        grad[xo] = Ms[m] * ga + Ms[m + 3] * gb + Ms[m + 6] * gc;
        grad[xo + 1] = Ms[m + 1] * ga + Ms[m + 4] * gb + Ms[m + 7] * gc;
        grad[xo + 2] = Ms[m + 2] * ga + Ms[m + 5] * gb + Ms[m + 8] * gc;
        // the rotation vector: Jl^T det(M) M^T (dE/dtheta)
        const ta = gU[off + o], tb = gU[off + o + 1], tc = gU[off + o + 2];
        const det = det9(Ms, m);
        const h0 = det * (Ms[m] * ta + Ms[m + 3] * tb + Ms[m + 6] * tc);
        const h1 = det * (Ms[m + 1] * ta + Ms[m + 4] * tb + Ms[m + 7] * tc);
        const h2 = det * (Ms[m + 2] * ta + Ms[m + 5] * tb + Ms[m + 8] * tc);
        leftJacobian(x[xo + 3], x[xo + 4], x[xo + 5], jac);
        grad[xo + 3] = jac[0] * h0 + jac[3] * h1 + jac[6] * h2;
        grad[xo + 4] = jac[1] * h0 + jac[4] * h1 + jac[7] * h2;
        grad[xo + 5] = jac[2] * h0 + jac[5] * h1 + jac[8] * h2;
        // the centre's motion with the size: dy/ds = (y - M v) / s (y / s for a piece at its anchor)
        const free = stuck[q] ? 0 : 1;
        const v0 = free && x[xo];
        const v1 = free && x[xo + 1];
        const v2 = free && x[xo + 2];
        const d0 = st.U[o] - (Ms[m] * v0 + Ms[m + 1] * v1 + Ms[m + 2] * v2);
        const d1 = st.U[o + 1] - (Ms[m + 3] * v0 + Ms[m + 4] * v1 + Ms[m + 5] * v2);
        const d2 = st.U[o + 2] - (Ms[m + 6] * v0 + Ms[m + 7] * v1 + Ms[m + 8] * v2);
        dsize += (ga * d0 + gb * d1 + gc * d2) / size;
      }
      grad[n] = dsize;
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
      const copy: State = { cells: Int32Array.from(st.cells), U: Float64Array.from(st.U), R: Float64Array.from(st.R!) };
      const f = finalizeState(model, copy, s.size, accept, { shrink: true });
      return f ? { config: f.config, score: f.score } : null;
    },

    validate(config): Validation {
      return serverRule(prob, config, serverBound);
    },

    seeds(rng: Rng, count: number) {
      const out: KernelState[] = [];
      for (let t = 0; t < count; t++) {
        // the cells themselves first (when the solid is the container's cell), then sequential random layouts
        const built = t === 0 ? model.constructedStart(rng) : null;
        const b = built ?? sequentialLayout(model, rng);
        out.push(stateOf(b.st, b.s));
      }
      return out;
    },

    swap(s, i, j) {
      if (i === j) return;
      const a = anchorsOf(s);
      const st = place(s);
      const Oi = Array.from(st.R!.subarray(9 * i, 9 * i + 9));
      const Oj = Array.from(st.R!.subarray(9 * j, 9 * j + 9));
      const wi = rotationVectorFor(Ms, 9 * i, a.R, 9 * i, Oj);
      const wj = rotationVectorFor(Ms, 9 * j, a.R, 9 * j, Oi);
      s.x.set(wi, dim * i + 3);
      s.x.set(wj, dim * j + 3);
    },
  };
  return kernel;
}

/**
 * A layout built piece by piece: each piece at a random cell, point and orientation, kept when the model finds no overlap of that piece
 * with those already placed or with an edge; the edge is raised by 6% after every failed round (up to 60 rounds, then the last, loosest
 * layout is returned and finish() judges it).
 */
export function sequentialLayout(model: SolidModel, rand: () => number): { st: State; s: number } {
  const poly = model.poly;
  const P = model.P;
  const C = poly.cells.length;
  const cells = new Int32Array(P);
  const U = new Float64Array(3 * P);
  const R = new Float64Array(9 * P);
  let s = model.startEdge();
  for (let round = 0; round < 60; round++, s *= 1.06) {
    let placedCount = 0;
    for (let tries = 0; placedCount < P && tries < 300 * P; tries++) {
      cells[placedCount] = randInt(rand, C);
      const p = randomCellPoint(poly, cells[placedCount], s, rand, 0.7);
      for (let k = 0; k < 3; k++) U[3 * placedCount + k] = p[k];
      R.set(randomRotation(rand), 9 * placedCount);
      const v = model.evaluate(cells.subarray(0, placedCount + 1), U.subarray(0, 3 * placedCount + 3), R.subarray(0, 9 * placedCount + 9), s, null, true, { focus: placedCount });
      if (v <= 0) placedCount++;
    }
    if (placedCount === P) return { st: { cells, U, R }, s };
  }
  return { st: { cells, U, R }, s };
}
