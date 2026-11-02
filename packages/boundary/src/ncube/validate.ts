/**
 * Float validator of the d-cube problems (d = 5..10, docs/spaces/BOUNDARY.md section 4B), tolerance tol (default 1e-9).
 *
 * The rule is the tesseract's (validate.ts, BOUNDARY.md 3): centres in their closed anchor cells (within tol), no hinge
 * (a (d-3)-face) inside a piece, no piece overlapping itself, interiors of distinct pieces disjoint. The d-cube's
 * lattice makes it exact in closed form (lattice.ts):
 *  - Balls (diameter 1): every centre at geodesic distance >= 1/2 - tol from the hinges (closed form), every two centres
 *    >= 1 - tol apart along every realised word, every geodesic loop at a centre >= 1 - tol. The words are found by reach
 *    (lattice.ts), for the pairs whose centres are < 1 apart in R^d only (the geodesic is never shorter than the chord).
 *  - Axis-aligned cubes (edge 1, no rotation): the cube centred at u is the box u +- 1/2. It crosses axis k when
 *    |u_k| + 1/2 > h + tol; valid alone exactly when it crosses at most one axis (4.3: two crossed axes put the anchor's
 *    own hinge inside it); its parts are the box clipped to the anchor and, beyond a crossed face, a slab in the
 *    neighbour; parts of distinct cubes in a cell must not overlap by more than tol in all axes. A rotated cube is not
 *    a d-cube problem here (`q` is refused).
 * No ridge snap is needed: every gluing is an exact signed permutation, so a piece is developed without rounding.
 *
 * Work is capped. A check stops after opts.maxWork search steps (default NCUBE_MAX_WORK) and refuses the configuration
 * with a plain reason, so a packing that crowds many balls against the hinges of a high-dimensional cube is refused
 * (sound, incomplete); the pair stage bins centres (balls) or the parts' corners (cubes) in a uniform grid, so the work
 * is linear in P, and it stops after maxReasons reasons.
 *
 * Untrusted input: numbers are read by index (holes are rejected) and every rejecting comparison is written so that NaN
 * rejects.
 */
import type { BoundaryVerdict } from "../types";
import { MAX_S } from "../validate";
import { unitBallVolume } from "./bounds";
import {
  NCubeName,
  WordLimitError,
  axesOf,
  embedCell,
  isNCubeName,
  ncubeClearance,
  ncubeDim,
  newBudget,
  realisedWords,
  startMap,
  stepMap,
} from "./lattice";

export interface NCubeValidateOptions {
  tol?: number;
  maxReasons?: number;
  /** Search steps before the check gives up (default NCUBE_MAX_WORK). */
  maxWork?: number;
}

/** Search steps (lattice.ts) one check may spend: about a second on one core ([CALC] BOUNDARY.md 4B). */
export const NCUBE_MAX_WORK = 4_000_000;

export interface NCubeProblem {
  polytope: NCubeName;
  piece: "ball" | "cube";
  P: number;
}

const fmt = (x: number) => (Number.isFinite(x) ? String(Number(x.toPrecision(10))) : String(x));

function isObject(x: unknown): x is Record<string, unknown> {
  return typeof x === "object" && x !== null && !Array.isArray(x);
}

function finiteNumbers(x: unknown, n: number): number[] | null {
  if (!Array.isArray(x) || x.length !== n) return null;
  const out: number[] = [];
  for (let k = 0; k < n; k++) {
    const v: unknown = x[k];
    if (typeof v !== "number" || !Number.isFinite(v)) return null;
    out.push(v);
  }
  return out;
}

export interface NCubeConfig {
  s: number;
  pieces: { cell: number; u: number[] }[];
}

/** Checks the shape of a configuration; returns the typed pieces or pushes reasons. */
export function readNCubeConfig(
  problem: NCubeProblem,
  config: unknown,
  tol: number,
  reasons: string[],
  maxReasons: number,
): NCubeConfig | null {
  const push = (r: string) => {
    if (reasons.length < maxReasons) reasons.push(r);
  };
  const d = ncubeDim(problem.polytope);
  if (!d) {
    push(`Unknown polytope "${String(problem.polytope)}".`);
    return null;
  }
  const n = d - 1;
  if (!isObject(config)) {
    push("The configuration must be an object with s and pieces.");
    return null;
  }
  const s = config.s;
  if (typeof s !== "number" || !Number.isFinite(s)) push("The edge s must be a finite number.");
  else if (!(s > 0 && s <= MAX_S)) push(`The edge s must be greater than 0 and at most ${MAX_S}; it is ${fmt(s)}.`);
  const list = config.pieces;
  if (!Array.isArray(list)) {
    push("The configuration needs a list of pieces.");
    return null;
  }
  if (list.length !== problem.P)
    push(`The configuration has ${list.length} pieces; the problem needs exactly ${problem.P}.`);
  if (reasons.length) return null;
  const h = (s as number) / 2;
  const pieces: NCubeConfig["pieces"] = [];
  for (let i = 0; i < list.length; i++) {
    const raw: unknown = list[i];
    if (!isObject(raw)) {
      push(`Piece ${i} must be an object with cell and u.`);
      continue;
    }
    let ok = true;
    const cell = raw.cell;
    if (!(typeof cell === "number" && Number.isInteger(cell) && cell >= 0 && cell < 2 * d)) {
      push(`Piece ${i}: its cell must be a whole number from 0 to ${2 * d - 1}.`);
      ok = false;
    }
    const u = finiteNumbers(raw.u, n);
    if (!u) {
      push(`Piece ${i}: its centre u must be ${n} finite numbers.`);
      ok = false;
    }
    if (raw.q !== undefined) {
      push(`Piece ${i}: a piece on the ${problem.polytope} has no rotation q (cubes stay axis-aligned).`);
      ok = false;
    }
    if (!ok || !u) continue;
    const ex = Math.max(...u.map((x) => Math.abs(x))) - h;
    if (!(ex <= tol)) push(`Piece ${i}: its centre lies outside its anchor cell ${cell} by ${fmt(ex)}.`);
    pieces.push({ cell: cell as number, u });
  }
  if (reasons.length) return null;
  return { s: s as number, pieces };
}

/** Pairs of pieces by index whose R^d centres are less than `reach` apart (a grid on the first four coordinates). */
export function nearPairs(
  d: number,
  cells: number[],
  U: Float64Array,
  n: number,
  h: number,
  reach: number,
): [number, number][] {
  const P = cells.length;
  const X = cells.map((c, i) => embedCell(d, c, U, i * n, h));
  const g = new Map<string, number[]>();
  const key = (x: Float64Array, a: number, b: number, c: number, e: number) =>
    `${Math.floor(x[0] / reach) + a},${Math.floor(x[1] / reach) + b},${Math.floor(x[2] / reach) + c},${Math.floor(x[3] / reach) + e}`;
  for (let i = 0; i < P; i++) {
    const k = key(X[i], 0, 0, 0, 0);
    const l = g.get(k);
    if (l) l.push(i);
    else g.set(k, [i]);
  }
  const out: [number, number][] = [];
  const r2 = reach * reach;
  for (let i = 0; i < P; i++)
    for (let a = -1; a <= 1; a++)
      for (let b = -1; b <= 1; b++)
        for (let c = -1; c <= 1; c++)
          for (let e = -1; e <= 1; e++)
            for (const j of g.get(key(X[i], a, b, c, e)) ?? []) {
              if (j <= i) continue;
              let q = 0;
              for (let t = 0; t < d; t++) q += (X[i][t] - X[j][t]) ** 2;
              if (q < r2) out.push([i, j]);
            }
  return out;
}

function validateBalls(
  d: number,
  name: string,
  cfg: NCubeConfig,
  tol: number,
  maxReasons: number,
  maxWork: number,
): string[] {
  const reasons: string[] = [];
  const push = (r: string) => {
    if (reasons.length < maxReasons) reasons.push(r);
  };
  const n = d - 1;
  const { s, pieces } = cfg;
  const h = s / 2;
  const P = pieces.length;
  const U = new Float64Array(n * P);
  const cells = pieces.map((p, i) => {
    U.set(p.u, i * n);
    return p.cell;
  });
  pieces.forEach((_, i) => {
    const c = ncubeClearance(U, i * n, n, h);
    if (!(c >= 0.5 - tol))
      push(
        `Ball ${i} is ${fmt(c)} from a (d-3)-face of the ${name}; it must be at least 0.5, or the face passes through the ball.`,
      );
  });
  if (reasons.length) return reasons;
  const budget = newBudget(maxWork);
  try {
    const limit = 1 - tol;
    for (const [i, j] of nearPairs(d, cells, U, n, h, 1)) {
      if (reasons.length >= maxReasons) break;
      let best = Infinity;
      for (const r of realisedWords(d, s, cells[i], U, i * n, cells[j], U, j * n, budget)) best = Math.min(best, r.len);
      if (!(best >= limit))
        push(
          `Balls ${i} and ${j} overlap: their centres are ${fmt(best)} apart along the boundary; they must be at least 1 apart.`,
        );
    }
    for (let i = 0; i < P && reasons.length < maxReasons; i++) {
      let best = Infinity;
      for (const r of realisedWords(d, s, cells[i], U, i * n, cells[i], U, i * n, budget, { excludeEmpty: true }))
        if (r.len > 1e-12) best = Math.min(best, r.len);
      if (!(best >= limit))
        push(
          `Ball ${i} overlaps itself: a geodesic loop of length ${fmt(best)} leaves and returns to its centre; it must be at least 1.`,
        );
    }
  } catch (e) {
    if (!(e instanceof WordLimitError)) throw e;
    push(
      `The check reached its work limit (${e.limit} search steps) before it could decide: too many balls crowd the (d-3)-faces of the ${name}.`,
    );
  }
  return reasons;
}

/** A part of a cube: a box in a cell's local coordinates. */
interface Part {
  piece: number;
  cell: number;
  lo: Float64Array;
  hi: Float64Array;
}

/**
 * The parts of the axis-aligned cube centred at u (anchor cell c) as boxes, and the axes it crosses (depth > tol).
 * With one crossed axis k (side tau) the slice beyond the face y_k = tau h, y_k in [h, u_k + 1/2] (tau > 0), has the
 * other axes' bounds, and the one-crossing word maps it to a box of the neighbour cell (j, tau'); the old normal axis i
 * becomes a tangent axis of the neighbour, of depth d_k into it.
 */
export function cubeParts(
  d: number,
  s: number,
  piece: number,
  c: number,
  u: ArrayLike<number>,
  o: number,
  tol: number,
): { parts: Part[]; crossed: number[]; depth: number[] } {
  const n = d - 1;
  const h = s / 2;
  const depth: number[] = [];
  const crossed: number[] = [];
  for (let k = 0; k < n; k++) {
    const dep = Math.abs(u[o + k]) + 0.5 - h;
    depth.push(dep);
    if (dep > tol) crossed.push(k);
  }
  const lo = new Float64Array(n);
  const hi = new Float64Array(n);
  for (let k = 0; k < n; k++) {
    lo[k] = u[o + k] - 0.5;
    hi[k] = u[o + k] + 0.5;
    if (depth[k] > tol) {
      if (u[o + k] > 0) hi[k] = h;
      else lo[k] = -h;
    }
  }
  const parts: Part[] = [{ piece, cell: c, lo, hi }];
  if (crossed.length === 1) {
    const k = crossed[0];
    const tau = u[o + k] > 0 ? 1 : -1;
    const blo = Float64Array.from(lo);
    const bhi = Float64Array.from(hi);
    if (tau > 0) {
      blo[k] = h;
      bhi[k] = u[o + k] + 0.5;
    } else {
      blo[k] = u[o + k] - 0.5;
      bhi[k] = -h;
    }
    const w = stepMap(d, startMap(d, c), k, tau);
    // x = A y + h b for the two corners, then the neighbour's local axes
    const ax2 = axesOf(d, w.cell >> 1);
    const A = (y: Float64Array) => {
      const x = new Float64Array(d);
      for (let r = 0; r < n; r++) x[w.rows[r]] += w.sgn[r] * y[r];
      for (let t = 0; t < d; t++) x[t] += h * w.b[t];
      return x;
    };
    const x0 = A(blo);
    const x1 = A(bhi);
    const nlo = new Float64Array(n);
    const nhi = new Float64Array(n);
    for (let q = 0; q < n; q++) {
      nlo[q] = Math.min(x0[ax2[q]], x1[ax2[q]]);
      nhi[q] = Math.max(x0[ax2[q]], x1[ax2[q]]);
    }
    parts.push({ piece, cell: w.cell, lo: nlo, hi: nhi });
  }
  return { parts, crossed, depth };
}

function validateCubes(d: number, name: string, cfg: NCubeConfig, tol: number, maxReasons: number): string[] {
  const reasons: string[] = [];
  const push = (r: string) => {
    if (reasons.length < maxReasons) reasons.push(r);
  };
  const n = d - 1;
  const { s, pieces } = cfg;
  const h = s / 2;
  const U = new Float64Array(n * pieces.length);
  pieces.forEach((p, i) => U.set(p.u, i * n));
  const byCell: Part[][] = Array.from({ length: 2 * d }, () => []);
  for (let i = 0; i < pieces.length; i++) {
    const { parts, crossed, depth } = cubeParts(d, s, i, pieces[i].cell, U, i * n, tol);
    if (crossed.length >= 2) {
      const worst = depth.slice().sort((a, b) => b - a)[1];
      push(
        `Cube ${i} crosses the faces of its cell in ${crossed.length} axes, so a (d-3)-face of the ${name} passes through its interior (depth ${fmt(worst)}).`,
      );
      continue;
    }
    for (const p of parts) byCell[p.cell].push(p);
  }
  if (reasons.length) return reasons;
  // the parts' low corners in a grid of side 2 on the first three axes: a part is at most 1 long on each axis
  const G = 2;
  const kc = (x: number) => Math.floor(x / G);
  void h;
  byCell.forEach((list, cell) => {
    if (reasons.length >= maxReasons) return;
    const grid = new Map<string, number[]>();
    const key = (a: number, b: number, c: number) => `${a},${b},${c}`;
    list.forEach((p, a) => {
      const k = key(kc(p.lo[0]), kc(p.lo[1]), kc(p.lo[2]));
      const l = grid.get(k);
      if (l) l.push(a);
      else grid.set(k, [a]);
    });
    for (let a = 0; a < list.length && reasons.length < maxReasons; a++) {
      const A = list[a];
      for (let dx = -1; dx <= 1; dx++)
        for (let dy = -1; dy <= 1; dy++)
          for (let dz = -1; dz <= 1; dz++)
            for (const b of grid.get(key(kc(A.lo[0]) + dx, kc(A.lo[1]) + dy, kc(A.lo[2]) + dz)) ?? []) {
              if (b <= a) continue;
              const B = list[b];
              // the least overlap over the axes is the penetration depth; all overlaps must exceed tol
              let depth = Infinity;
              for (let k = 0; k < n; k++) {
                const ov = Math.min(A.hi[k], B.hi[k]) - Math.max(A.lo[k], B.lo[k]);
                if (ov < depth) depth = ov;
                if (!(depth > tol)) break;
              }
              if (!(depth > tol)) continue;
              const [i, j] = A.piece < B.piece ? [A.piece, B.piece] : [B.piece, A.piece];
              push(`Cubes ${i} and ${j} overlap in cell ${cell} (depth ${fmt(depth)}).`);
            }
    }
  });
  return reasons;
}

/**
 * Validates a configuration (untrusted input) for a d-cube problem: ok with the score s, or the reasons as plain
 * sentences. Throws a RangeError on a bad tol (a caller's mistake).
 */
export function validateNCube(
  problem: NCubeProblem,
  config: unknown,
  opts: NCubeValidateOptions = {},
): BoundaryVerdict {
  const tol = opts.tol ?? 1e-9;
  if (!(typeof tol === "number" && Number.isFinite(tol) && tol >= 0))
    throw new RangeError(`validateNCube: tol must be a finite number >= 0, not ${String(tol)}`);
  const mr = opts.maxReasons ?? 20;
  const maxReasons = mr >= 1 ? Math.floor(mr) : 1;
  const mw = opts.maxWork ?? NCUBE_MAX_WORK;
  const maxWork = mw >= 1 ? Math.floor(mw) : 1;
  const reasons: string[] = [];
  const bad = (): BoundaryVerdict => ({
    ok: false,
    score: null,
    reasons: reasons.length ? reasons : ["The configuration is malformed."],
  });
  if (!isObject(problem) || !isNCubeName(problem.polytope)) {
    reasons.push(`Unknown polytope "${String((problem as { polytope?: unknown })?.polytope)}".`);
    return bad();
  }
  if (problem.piece !== "ball" && problem.piece !== "cube") {
    reasons.push(`Unknown piece "${String(problem.piece)}".`);
    return bad();
  }
  if (!(typeof problem.P === "number" && Number.isInteger(problem.P) && problem.P >= 1)) {
    reasons.push("The problem's P must be a whole number from 1.");
    return bad();
  }
  const d = ncubeDim(problem.polytope);
  const cfg = readNCubeConfig(problem, config, tol, reasons, maxReasons);
  if (!cfg) return bad();
  const n = d - 1;
  // volume bound [CALC]: the pieces' volumes add up to at most the boundary's, 2d s^n (1e-6 absorbs overlaps within tol)
  const pieceVolume = problem.piece === "ball" ? unitBallVolume(n) / 2 ** n : 1;
  const volume = 2 * d * cfg.s ** n;
  if (!(cfg.pieces.length * pieceVolume <= volume * (1 + 1e-6)))
    return {
      ok: false,
      score: null,
      reasons: [
        `${cfg.pieces.length} ${problem.piece}s have a total volume of ${fmt(cfg.pieces.length * pieceVolume)}, more than the ${fmt(volume)} of the ${problem.polytope}'s boundary at s = ${fmt(cfg.s)}.`,
      ],
    };
  const geo =
    problem.piece === "ball"
      ? validateBalls(d, problem.polytope, cfg, tol, maxReasons, maxWork)
      : validateCubes(d, problem.polytope, cfg, tol, maxReasons);
  return geo.length ? { ok: false, score: null, reasons: geo } : { ok: true, score: cfg.s, reasons: [] };
}
