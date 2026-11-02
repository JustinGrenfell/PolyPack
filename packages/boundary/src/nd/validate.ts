/**
 * Float validator of the simplex, orthoplex and demicube problems (d = 5..10, balls of diameter 1; docs/spaces/BOUNDARY.md section 4C),
 * tolerance tol (default 1e-9). The rule is 4.2's, which is dimension-free: a configuration is valid exactly when
 *  - every centre lies in its closed anchor cell (within tol),
 *  - every centre is at geodesic distance >= 1/2 - tol from the hinges (the (d-3)-faces),
 *  - every two centres are >= 1 - tol apart along the boundary,
 *  - no geodesic loop at a centre is shorter than 1 - tol.
 * The engine is chains.ts: half-chains of reach 1/2 from every centre, hinge distances inside their windows, and joins of two
 * half-chains in a common cell for pairs and loops. Pairs are searched only when the centres are < 1 apart in R^D (a geodesic is never
 * shorter than the chord).
 *
 * A centre within the snap distance eta of a ridge plane is moved onto it (the move is charged to the verdict: thresholds of that piece
 * are tightened by the distance it moved) and then belongs to both cells, which are both used as anchors, as in the 4-polytope engine.
 *
 * Work is capped. A check stops after opts.maxWork steps (default ND_MAX_WORK) and refuses the configuration with a plain reason: a
 * packing that crowds balls against the hinges, or a container too small for its balls, is refused (sound, incomplete), never accepted.
 *
 * Untrusted input: numbers are read by index (holes are rejected) and every rejecting comparison is written so that NaN rejects.
 */
import type { BoundaryVerdict } from "../types";
import { MAX_S } from "../validate";
import { type NDBudget, type NDNode, NDLimitError, develop, hingeDistance, joinDistance, newNDBudget } from "./chains";
import { type NDContainer, getND } from "./container";
import { type NDName, isNDName } from "./names";
import { ndBoundaryVolume, ndBallVolume } from "./bounds";

export interface NDValidateOptions {
  tol?: number;
  maxReasons?: number;
  /** Search steps before the check gives up (default ND_MAX_WORK). */
  maxWork?: number;
  /** Receives the steps used. */
  stats?: { steps: number };
}

/** Search steps one check may spend (chains.ts): BOUNDARY.md 4C.7 gives the measurements it comes from. */
export const ND_MAX_WORK = 1_500_000;

export interface NDProblem {
  polytope: NDName;
  piece: "ball";
  P: number;
}

export interface NDConfig {
  s: number;
  pieces: { cell: number; u: number[] }[];
}

const fmt = (x: number) => (Number.isFinite(x) ? String(Number(x.toPrecision(10))) : String(x));
const isObject = (x: unknown): x is Record<string, unknown> => typeof x === "object" && x !== null && !Array.isArray(x);

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

/** Checks the shape of a configuration; returns the typed pieces or pushes reasons. */
export function readNDConfig(problem: NDProblem, config: unknown, tol: number, reasons: string[], maxReasons: number): NDConfig | null {
  const push = (r: string) => {
    if (reasons.length < maxReasons) reasons.push(r);
  };
  const C = getND(problem.polytope);
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
  if (list.length !== problem.P) push(`The configuration has ${list.length} pieces; the problem needs exactly ${problem.P}.`);
  if (reasons.length) return null;
  const pieces: NDConfig["pieces"] = [];
  for (let i = 0; i < list.length; i++) {
    const raw: unknown = list[i];
    if (!isObject(raw)) {
      push(`Piece ${i} must be an object with cell and u.`);
      continue;
    }
    let ok = true;
    const cell = raw.cell;
    if (!(typeof cell === "number" && Number.isInteger(cell) && cell >= 0 && cell < C.C)) {
      push(`Piece ${i}: its cell must be a whole number from 0 to ${C.C - 1}.`);
      ok = false;
    }
    const u = finiteNumbers(raw.u, C.n);
    if (!u) {
      push(`Piece ${i}: its centre u must be ${C.n} finite numbers.`);
      ok = false;
    }
    if (raw.q !== undefined) {
      push(`Piece ${i}: a ball on the ${problem.polytope} has no rotation q.`);
      ok = false;
    }
    if (!ok || !u) continue;
    pieces.push({ cell: cell as number, u });
  }
  if (reasons.length) return null;
  void tol;
  return { s: s as number, pieces };
}

/** The snap distance of a centre to a ridge plane (4-polytope engine's ridgeSnap, chains.ts there). */
export function ndSnap(s: number, tol: number): number {
  const base = 1e-14 * Math.max(1, s);
  return Math.max(2 ** -52 * Math.max(1, s), Math.min(base, tol / 10));
}

export interface NDAnchored {
  /** The claimed cell and the centre in ambient coordinates, clamped into it and snapped to the ridges within eta. */
  x: Float64Array;
  cell: number;
  /** How far the clamp and snap moved the centre. */
  moved: number;
  /** How far the centre was outside its closed cell before (<= 0 inside). */
  outside: number;
  /** All the cells the centre belongs to (the claimed one first) and the centre in each one's hyperplane. */
  anchors: { cell: number; x: Float64Array }[];
}

/** Reads a centre of cell c (local u at edge s): clamped into the closed cell, snapped onto the ridge planes within eta. */
export function anchorPiece(C: NDContainer, s: number, c: number, u: ArrayLike<number>, eta: number, uo = 0): NDAnchored {
  const D = C.D;
  const cl = C.cell(c);
  const x0 = C.toAmbient(c, u, s, uo);
  const x = Float64Array.from(x0);
  let outside = -Infinity;
  const F = cl.nbr.length;
  for (let k = 0; k < F; k++) {
    let v = -s * cl.hn[k];
    for (let j = 0; j < D; j++) v += cl.nn[k * D + j] * x[j];
    if (v > outside) outside = v;
  }
  const nu = new Float64Array(D);
  for (let pass = 0; pass < 6; pass++) {
    let moved = false;
    for (let k = 0; k < F; k++) {
      let v = -s * cl.hn[k];
      let c1 = 0;
      for (let j = 0; j < D; j++) {
        v += cl.nn[k * D + j] * x[j];
        c1 += cl.nn[k * D + j] * cl.n[j];
      }
      if (v > 0 || (v > -eta && v !== 0)) {
        const f = v / cl.nu2[k];
        for (let j = 0; j < D; j++) {
          nu[j] = cl.nn[k * D + j] - c1 * cl.n[j];
          x[j] -= f * nu[j];
        }
        moved = true;
      }
    }
    if (!moved) break;
  }
  let mv = 0;
  for (let j = 0; j < D; j++) mv += (x[j] - x0[j]) ** 2;
  const anchors: NDAnchored["anchors"] = [{ cell: c, x }];
  for (let k = 0; k < F; k++) {
    let v = -s * cl.hn[k];
    for (let j = 0; j < D; j++) v += cl.nn[k * D + j] * x[j];
    if (Math.abs(v) <= 4 * eta) {
      const g = C.cell(cl.nbr[k]);
      const y = Float64Array.from(x);
      let e = -s * g.h;
      for (let j = 0; j < D; j++) e += g.n[j] * y[j];
      for (let j = 0; j < D; j++) y[j] -= e * g.n[j];
      anchors.push({ cell: cl.nbr[k], x: y });
    }
  }
  return { x, cell: c, moved: Math.sqrt(mv), outside, anchors };
}

/**
 * Reach of the half-chains: half the diameter of a ball, with a hair for rounding. It must stay below the margins of the packings the server is to
 * accept (records carry margins of 1e-12 and more): a hinge or a pair within the reach makes a tie, and ties are where the search cost explodes (4C.8).
 */
export const ND_REACH = 0.5 * (1 + 1e-12);

function validateBalls(C: NDContainer, cfg: NDConfig, tol: number, maxReasons: number, budget: NDBudget): string[] {
  const reasons: string[] = [];
  const push = (r: string) => {
    if (reasons.length < maxReasons) reasons.push(r);
  };
  const { s, pieces } = cfg;
  const P = pieces.length;
  const eta = ndSnap(s, tol);
  const at: NDAnchored[] = pieces.map((p) => anchorPiece(C, s, p.cell, p.u, eta));
  const bad = new Uint8Array(P);
  at.forEach((a, i) => {
    if (!(a.outside <= tol)) {
      push(`Piece ${i}: its centre lies outside its anchor cell ${a.cell} by ${fmt(a.outside)}.`);
      bad[i] = 1;
    }
  });
  if (reasons.length) return reasons;
  // the half-chains of every centre, from every cell it belongs to
  const chains: NDNode[][][] = [];
  const hinge: number[] = [];
  for (let i = 0; i < P; i++) {
    // the centre's own cell first: a centre nearer than 1/2 to a hinge of it is refused before any chain is built
    const rootOnly = at[i].anchors.map((an) => develop(C, s, an.cell, an.x, 0, budget));
    let h0 = Infinity;
    for (const r of rootOnly) h0 = Math.min(h0, hingeDistance(C, s, r, ND_REACH, budget).dist);
    if (!(h0 >= 0.5 - tol + at[i].moved)) {
      push(`Ball ${i} is ${fmt(h0)} from a (d-3)-face of the ${C.name}; it must be at least 0.5, or the face passes through the ball.`);
      bad[i] = 1;
      chains.push([]);
      hinge.push(h0);
      continue;
    }
    const cs = at[i].anchors.map((an) => develop(C, s, an.cell, an.x, ND_REACH, budget));
    chains.push(cs);
    let h = Infinity;
    for (const c of cs) h = Math.min(h, hingeDistance(C, s, c, ND_REACH, budget).dist);
    hinge.push(h);
    if (!(h >= 0.5 - tol + at[i].moved))
      push(`Ball ${i} is ${fmt(h)} from a (d-3)-face of the ${C.name} along a path over a ridge; it must be at least 0.5, or the face passes through the ball.`);
  }
  if (reasons.length) return reasons;
  // pairs
  for (let i = 0; i < P && reasons.length < maxReasons; i++)
    for (let j = i + 1; j < P && reasons.length < maxReasons; j++) {
      let q = 0;
      for (let r = 0; r < C.D; r++) q += (at[i].x[r] - at[j].x[r]) ** 2;
      const lim = 1 - tol + at[i].moved + at[j].moved;
      if (!(q < 1 * 1)) continue;
      let best = Infinity;
      for (const A of chains[i]) for (const B of chains[j]) best = Math.min(best, joinDistance(C, s, A, B, lim, budget, false).dist);
      if (best < lim) push(`Balls ${i} and ${j} overlap: their centres are ${fmt(best)} apart along the boundary; they must be at least 1 apart.`);
    }
  // loops
  for (let i = 0; i < P && reasons.length < maxReasons; i++) {
    const lim = 1 - tol + 2 * at[i].moved;
    let best = Infinity;
    for (let a = 0; a < chains[i].length; a++)
      for (let b = a; b < chains[i].length; b++) best = Math.min(best, joinDistance(C, s, chains[i][a], chains[i][b], lim, budget, a === b).dist);
    if (best < lim) push(`Ball ${i} overlaps itself: a geodesic loop of length ${fmt(best)} leaves and returns to its centre; it must be at least 1.`);
  }
  return reasons;
}

/** Validates a configuration of balls on a simplex, orthoplex or demicube. Never throws on input. */
export function validateND(problem: NDProblem, config: unknown, opts: NDValidateOptions = {}): BoundaryVerdict {
  const tol = typeof opts.tol === "number" && opts.tol >= 0 ? opts.tol : 1e-9;
  const mr = opts.maxReasons ?? 20;
  const maxReasons = mr >= 1 ? Math.floor(mr) : 1;
  const mw = opts.maxWork ?? ND_MAX_WORK;
  const maxWork = mw >= 1 ? Math.floor(mw) : 1;
  const reasons: string[] = [];
  const bad = (): BoundaryVerdict => ({ ok: false, score: null, reasons: reasons.length ? reasons : ["The configuration is malformed."] });
  if (!isObject(problem) || !isNDName(problem.polytope)) {
    reasons.push(`Unknown polytope "${String((problem as { polytope?: unknown })?.polytope)}".`);
    return bad();
  }
  if (problem.piece !== "ball") {
    reasons.push(`Unknown piece "${String(problem.piece)}": the ${problem.polytope} takes balls only.`);
    return bad();
  }
  if (!(typeof problem.P === "number" && Number.isInteger(problem.P) && problem.P >= 1)) {
    reasons.push("The problem's P must be a whole number from 1.");
    return bad();
  }
  const C = getND(problem.polytope);
  const cfg = readNDConfig(problem, config, tol, reasons, maxReasons);
  if (!cfg) return bad();
  // volume bound [CALC]: the balls' volumes add up to at most the boundary's (1e-6 absorbs overlaps within tol)
  const pieceVolume = ndBallVolume(C.n);
  const volume = ndBoundaryVolume(C.family, C.d) * cfg.s ** C.n;
  if (!(cfg.pieces.length * pieceVolume <= volume * (1 + 1e-6)))
    return {
      ok: false,
      score: null,
      reasons: [
        `${cfg.pieces.length} balls have a total volume of ${fmt(cfg.pieces.length * pieceVolume)}, more than the ${fmt(volume)} of the ${problem.polytope}'s boundary at s = ${fmt(cfg.s)}.`,
      ],
    };
  const budget = newNDBudget(maxWork);
  let geo: string[];
  try {
    geo = validateBalls(C, cfg, tol, maxReasons, budget);
  } catch (e) {
    if (!(e instanceof NDLimitError)) throw e;
    geo = [`The check reached its work limit (${maxWork} search steps) before it could decide: too many balls crowd the (d-3)-faces of the ${problem.polytope}.`];
  }
  if (opts.stats) opts.stats.steps = budget.steps;
  return geo.length ? { ok: false, score: null, reasons: geo } : { ok: true, score: cfg.s, reasons: [] };
}

export interface NDMeasure {
  /** Geodesic distance of each centre to the hinges (Infinity when above the reach 1/2), and the 'crowded' balls (within 1/2 of a hinge of their own cell, not developed). */
  clearance: number[];
  /** The shortest straight path between each pair of centres closer than 1 in R^D (Infinity when none is shorter than 1), and the shortest loop at each centre. */
  pairs: { i: number; j: number; dist: number }[];
  loops: number[];
  steps: number;
}

/**
 * The distances the rule is made of, for the cross-check against the exact certifier and the docs' measurements (no verdict): the
 * clearance of each centre, the pair distances below 1 and the loops below 1. A centre within 1/2 of a hinge of its own cell is not
 * developed (its pairs and loop are reported as the chord and Infinity).
 */
export function measureND(problem: NDProblem, config: NDConfig, opts: { maxWork?: number } = {}): NDMeasure {
  const C = getND(problem.polytope);
  const { s, pieces } = config;
  const P = pieces.length;
  const budget = newNDBudget(opts.maxWork ?? ND_MAX_WORK);
  const eta = ndSnap(s, 0);
  const at = pieces.map((p) => anchorPiece(C, s, p.cell, p.u, eta));
  const chains: NDNode[][][] = [];
  const clearance: number[] = [];
  const crowded: boolean[] = [];
  for (let i = 0; i < P; i++) {
    const roots = at[i].anchors.map((an) => develop(C, s, an.cell, an.x, 0, budget));
    let h = Infinity;
    for (const r of roots) h = Math.min(h, hingeDistance(C, s, r, ND_REACH, budget).dist);
    if (h < 0.5) {
      crowded.push(true);
      chains.push([]);
      clearance.push(h);
      continue;
    }
    crowded.push(false);
    const cs = at[i].anchors.map((an) => develop(C, s, an.cell, an.x, ND_REACH, budget));
    chains.push(cs);
    for (const c of cs) h = Math.min(h, hingeDistance(C, s, c, ND_REACH, budget).dist);
    clearance.push(h);
  }
  const pairs: NDMeasure["pairs"] = [];
  for (let i = 0; i < P; i++)
    for (let j = i + 1; j < P; j++) {
      let q = 0;
      for (let r = 0; r < C.D; r++) q += (at[i].x[r] - at[j].x[r]) ** 2;
      if (!(q < 1)) continue;
      if (crowded[i] || crowded[j]) {
        pairs.push({ i, j, dist: Math.sqrt(q) });
        continue;
      }
      let best = Infinity;
      for (const A of chains[i]) for (const B of chains[j]) best = Math.min(best, joinDistance(C, s, A, B, Math.min(best, 1), budget, false).dist);
      pairs.push({ i, j, dist: best });
    }
  const loops: number[] = [];
  for (let i = 0; i < P; i++) {
    let best = Infinity;
    for (let a = 0; a < chains[i].length; a++)
      for (let b = a; b < chains[i].length; b++) best = Math.min(best, joinDistance(C, s, chains[i][a], chains[i][b], Math.min(best, 1), budget, a === b).dist);
    loops.push(best);
  }
  return { clearance, pairs, loops, steps: budget.steps };
}
