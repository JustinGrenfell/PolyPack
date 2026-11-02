/**
 * The server's rule for a d-cube configuration, for the generic methods' kernel (kernel.ts) and the spaces registry
 * (packages/spaces/src/validate.ts). In order, as the tesseract's (kernel-validate.ts):
 *  1. P is at most the problem set's largest, NCUBE_MAX_P;
 *  2. the form, in the server's words: s a finite number in (0, 1e4]; a list of P pieces {cell, u}, cell an integer from 0 to
 *     2d - 1, u of d - 1 finite numbers of at most 1e4 in absolute value; no rotation q (balls and cubes alike: cubes stay
 *     axis-aligned);
 *  3. s at least the problem's lower bound less 1e-9 relative;
 *  4. validateNCube at tolerance 1e-9 with the check's work cap and 20 reasons.
 */
import type { Validation } from "@polypack/methods";
import { MAX_S } from "../validate";
import { NCUBE_MAX_P, ncubeBound } from "./bounds";
import { NCubeName, ncubeDim } from "./lattice";
import { NCUBE_MAX_WORK, NCubeProblem, validateNCube } from "./validate";

export const NCUBE_SERVER_RULE = {
  tol: 1e-9,
  maxReasons: 20,
  maxCoord: 1e4,
  boundSlack: 1e-9,
  maxWork: NCUBE_MAX_WORK,
} as const;

const isObject = (v: unknown): v is Record<string, unknown> => typeof v === "object" && v !== null && !Array.isArray(v);

/** What is wrong with the form of a d-cube configuration, in the server's words; empty when it is well formed. */
export function ncubeFormErrors(problem: NCubeProblem, raw: unknown): string[] {
  const R = NCUBE_SERVER_RULE;
  if (!isObject(raw)) return ["config must be an object"];
  const d = ncubeDim(problem.polytope);
  const n = d - 1;
  const errors: string[] = [];
  const num = (v: unknown, path: string, max: number, positive = false): number => {
    if (typeof v !== "number" || !Number.isFinite(v)) {
      errors.push(`${path} must be a finite number`);
      return NaN;
    }
    if (Math.abs(v) > max) errors.push(`${path} must be at most ${max} in absolute value`);
    else if (positive && !(v > 0)) errors.push(`${path} must be positive`);
    return v;
  };
  num(raw.s, "s", MAX_S, true);
  const list = raw.pieces;
  if (!Array.isArray(list)) return [...errors, "pieces must be an array"];
  if (list.length !== problem.P)
    return [...errors, `expected ${problem.P} pieces for this problem, got ${list.length}`];
  for (let i = 0; i < list.length; i++) {
    if (errors.length >= R.maxReasons) break;
    if (!(i in list)) continue;
    const p: unknown = list[i];
    if (!isObject(p)) {
      errors.push(`pieces[${i}] must be an object {cell, u}`);
      continue;
    }
    if (!Number.isInteger(p.cell) || (p.cell as number) < 0 || (p.cell as number) >= 2 * d)
      errors.push(`pieces[${i}].cell must be an integer from 0 to ${2 * d - 1}`);
    if (!Array.isArray(p.u) || p.u.length !== n)
      errors.push(`pieces[${i}].u must be [u0, ..., u${n - 1}] (${n} numbers)`);
    else for (let k = 0; k < n; k++) num((p.u as unknown[])[k], `pieces[${i}].u[${k}]`, R.maxCoord);
    if (p.q !== undefined)
      errors.push(`pieces[${i}] has no rotation q (pieces on the ${problem.polytope} stay axis-aligned)`);
  }
  return errors;
}

/** The server's verdict on a configuration of the problem (header). */
export function ncubeServerRule(problem: NCubeProblem, config: unknown, bound: number): Validation {
  const R = NCUBE_SERVER_RULE;
  const fail = (reasons: string[], score: number | null = null): Validation => ({ valid: false, score, reasons });
  const max = NCUBE_MAX_P[problem.polytope as NCubeName];
  if (!(problem.P <= max))
    return fail([`P = ${problem.P} is above the largest P of the ${problem.polytope} problems (${max})`]);
  const errors = ncubeFormErrors(problem, config);
  if (errors.length) return fail(errors.slice(0, R.maxReasons));
  const s = (config as { s: number }).s;
  if (s < bound * (1 - R.boundSlack))
    return fail([`the size ${s} is below the lower bound ${bound} of this problem`], s);
  let reasons: string[];
  try {
    const v = validateNCube(problem, config, { tol: R.tol, maxReasons: R.maxReasons, maxWork: R.maxWork });
    reasons = v.ok ? [] : v.reasons.length ? v.reasons.slice(0, R.maxReasons) : ["not a valid packing"];
  } catch (e) {
    reasons = [`the check failed: ${e instanceof Error ? e.message : String(e)}`];
  }
  return reasons.length ? fail(reasons, s) : { valid: true, score: s, reasons: [] };
}

/** The bound of a problem. */
export function ncubeProblemBoundOf(problem: NCubeProblem): number {
  return ncubeBound(ncubeDim(problem.polytope), problem.piece, problem.P);
}
