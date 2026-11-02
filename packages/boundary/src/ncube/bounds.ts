/**
 * Lower bounds, sizes and costs of the d-cube problems (docs/spaces/BOUNDARY.md section 4B).
 *
 * The container is the boundary of the d-cube [-s/2, s/2]^d: 2d facets, each a (d-1)-cube of edge s, so its (d-1)-volume
 * is 2d s^(d-1). Balls have diameter 1 and cubes edge 1. Both bounds below are proved and stored as the largest binary64
 * value not above the exact one (decided in integer arithmetic), as the tesseract's (kernel-validate.ts boundaryBound):
 *  - volume: the developed pieces have disjoint interiors and a development is injective on interiors, so
 *    P vol(piece) <= 2d s^(d-1). A ball's volume is V_n / 2^n, V_n the volume of the unit n-ball, n = d - 1;
 *  - one piece: a ball's centre is at distance >= 1/2 from the hinges, and the clearance is at most h sqrt 2 (the two
 *    slacks h - |u_k| are at most h), so s >= 1/sqrt 2; a cube needs s >= 1 (with s < 1 it crosses every axis of its cell,
 *    and two crossed axes put a hinge inside it, 4.3).
 */
import { NCUBE_NAMES, ncubeDim } from "./lattice";

/** Largest P of each d-cube problem set: the first tilings (P = 2d k^(d-1) at s = k) for d = 5, 6, then a common cap. */
export const NCUBE_MAX_P: Record<(typeof NCUBE_NAMES)[number], number> = {
  "5cube": 160,
  "6cube": 384,
  "7cube": 400,
  "8cube": 400,
  "9cube": 400,
  "10cube": 400,
};

/** floor(pi 10^60). */
const PI_60 = 3141592653589793238462643383279502884197169399375105820974944n;
const TEN60 = 10n ** 60n;

/** V_n = (num / den) pi^k for n = 4..9 (the volume of the unit n-ball). */
const BALL_VOLUME: Record<number, { num: bigint; den: bigint; k: number }> = {
  4: { num: 1n, den: 2n, k: 2 },
  5: { num: 8n, den: 15n, k: 2 },
  6: { num: 1n, den: 6n, k: 3 },
  7: { num: 16n, den: 105n, k: 3 },
  8: { num: 1n, den: 24n, k: 4 },
  9: { num: 32n, den: 945n, k: 4 },
};

/** The sign of x^p - num / den, exactly, for a positive finite x and positive integers num, den. */
function compare(x: number, p: number, num: bigint, den: bigint): number {
  const b = new DataView(new ArrayBuffer(8));
  b.setFloat64(0, x);
  const bits = b.getBigUint64(0);
  const exp = Number((bits >> 52n) & 0x7ffn);
  const frac = bits & ((1n << 52n) - 1n);
  const m = exp === 0 ? frac : frac | (1n << 52n);
  const e = (exp === 0 ? 1 : exp) - 1075;
  let left = m ** BigInt(p) * den;
  let right = num;
  if (e * p >= 0) left <<= BigInt(e * p);
  else right <<= BigInt(-e * p);
  return left > right ? 1 : left < right ? -1 : 0;
}

function step(x: number, dir: 1 | -1): number {
  const b = new DataView(new ArrayBuffer(8));
  b.setFloat64(0, x);
  b.setBigUint64(0, b.getBigUint64(0) + BigInt(dir));
  return b.getFloat64(0);
}

function largestFloat(check: (x: number) => boolean, approx: number): number {
  let f = approx;
  while (!check(f)) f = step(f, -1);
  while (check(step(f, 1))) f = step(f, 1);
  return f;
}

/** The bound on s of P pieces on the boundary of the d-cube (header). */
export function ncubeBound(d: number, piece: "ball" | "cube", P: number): number {
  const n = d - 1;
  const one = piece === "ball" ? largestFloat((x) => compare(x, 2, 1n, 2n) <= 0, Math.SQRT1_2) : 1;
  let volume: number;
  if (piece === "cube") {
    // P <= 2d s^n
    volume = largestFloat((x) => compare(x, n, BigInt(P), BigInt(2 * d)) <= 0, Math.pow(P / (2 * d), 1 / n));
  } else {
    // P V_n 2^-n <= 2d s^n with pi replaced by its truncation (a smaller bound, so a valid one)
    const v = BALL_VOLUME[n];
    const num = BigInt(P) * v.num * PI_60 ** BigInt(v.k);
    const den = v.den * 2n ** BigInt(n + 1) * BigInt(d) * TEN60 ** BigInt(v.k);
    volume = largestFloat(
      (x) => compare(x, n, num, den) <= 0,
      Math.pow((P * Math.pow(Math.PI, v.k) * Number(v.num)) / (Number(v.den) * 2 ** (n + 1) * d), 1 / n),
    );
  }
  return Math.max(one, volume);
}

/** The volume of the unit n-ball (floats), for the validator's volume test. */
export function unitBallVolume(n: number): number {
  const v = BALL_VOLUME[n];
  return (Number(v.num) / Number(v.den)) * Math.pow(Math.PI, v.k);
}

/** The bound of a problem by its container name. */
export function ncubeProblemBound(polytope: string, piece: "ball" | "cube", P: number): number {
  const d = ncubeDim(polytope);
  if (!d) throw new Error(`ncubeProblemBound: ${polytope} is no d-cube`);
  return ncubeBound(d, piece, P);
}

/** The titles of the d-cube containers. */
export const NCUBE_TITLES = {
  "5cube": "5-cube",
  "6cube": "6-cube",
  "7cube": "7-cube",
  "8cube": "8-cube",
  "9cube": "9-cube",
  "10cube": "10-cube",
} as const;

/**
 * Estimated milliseconds of the server's check per piece ([CALC] BOUNDARY.md 4B, 2026-11-01, one core): the cost of the
 * pair stage per ball measured on the cell layouts (one to twenty balls per cell, tight at their edges) was at most 0.5 ms
 * for the smallest P of d = 10 and 0.03 to 0.05 ms for the largest; the table is that worst per-piece cost times 3 for
 * crowded layouts and a cold engine. Cubes cost 0.02 ms each (boxes in a grid). The work of a check is also capped
 * (NCUBE_MAX_WORK steps, about a second), so a crowded layout is refused, not slow.
 */
export const NCUBE_BALL_MS: Record<(typeof NCUBE_NAMES)[number], number> = {
  "5cube": 0.3,
  "6cube": 0.35,
  "7cube": 0.5,
  "8cube": 0.6,
  "9cube": 1,
  "10cube": 1.5,
};

export function ncubeCheckMs(polytope: (typeof NCUBE_NAMES)[number], piece: "ball" | "cube", P: number): number {
  return piece === "cube" ? 0.02 * P : NCUBE_BALL_MS[polytope] * P;
}
