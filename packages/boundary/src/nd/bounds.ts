/**
 * Lower bounds, volumes and costs of the simplex, orthoplex and demicube problems (docs/spaces/BOUNDARY.md section 4C.5).
 *
 * Both bounds are proved and stored as the largest binary64 value not above the exact one, decided in integer arithmetic (as the d-cubes'
 * ncubeBound): the larger of
 *  - the volume bound P V_n 2^-n <= s^n vol_1: the developed balls have disjoint interiors, each is a Euclidean ball of radius 1/2 (its
 *    centre is 1/2 from the hinges, so the ball is the flat image of the Euclidean one), and the developments are injective on interiors,
 *    so their volumes add up to at most the boundary's; vol_1 is the boundary's n-volume at unit edge, n = d - 1 (formulas below);
 *  - the one-ball bound s >= 1 / (2 rho): a centre needs geodesic distance >= 1/2 from the hinges, hence from the hinges of its own
 *    cell, and every point of a cell lies within rho of a hinge of the cell, rho = the largest distance from the cell's centre O to a hinge
 *    of the cell. Proof (4A.4, for cells whose hinges are not all equidistant from O): a point p of the cell lies on a segment from O to a
 *    point q of a facet F of the cell; q is within the inradius r_F of the boundary of F (the farthest point of a convex polytope with
 *    a centre of symmetry group fixing it from its boundary is the centre), so within r_F of a facet e of F, a hinge of the cell; the distance
 *    to the convex set e is convex, so d(p, e) <= (1 - t) d(O, e) + t d(q, e) <= max(d(O, e), r_F) <= rho, as r_F <= d(O_F, e) <= d(O, e)
 *    for the centre O_F of F (the foot of O on F's plane). So s rho >= 1/2.
 *      simplex, orthoplex: cells are regular n-simplices of edge s, hinge distance s / sqrt(n^2 - 1) = s / sqrt(d (d - 2)): s >= sqrt(d (d-2)) / 2;
 *      demicube: the cells are regular n-simplices (rho^2 = 1/(n^2 - 1)) and n-demicubes of edge s, whose hinges between a demicube facet and a
 *      simplex facet are at squared distance (n^2 - 5 n + 8) / (8 (n - 1)) s^2 from the centre (and 1/4 s^2 between two demicube facets,
 *      n >= 5; not larger), so s^2 >= 2 (n - 1) / (n^2 - 5 n + 8).
 *
 * Volumes at unit edge, n = d - 1: a regular n-simplex sqrt(n+1) / (n! 2^(n/2)); an n-demicube 2^(-n/2) (1 - 2^(n-1)/n!) (the n-cube of
 * edge 1/sqrt 2 less 2^(n-1) corner simplices, each of volume (1/sqrt 2)^n / n!). The simplex's boundary has d + 1 simplices, the orthoplex's
 * 2^d, the demicube's 2^(d-1) simplices and 2d demicubes.
 */
import { type NDFamily, ND_NAMES, type NDName, ndParse } from "./names";

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

const factorial = (n: number): bigint => {
  let f = 1n;
  for (let i = 2; i <= n; i++) f *= BigInt(i);
  return f;
};

/** floor(sqrt(x)) of a non-negative BigInt (Newton's iteration). */
function isqrt(x: bigint): bigint {
  if (x < 2n) return x;
  let r = 1n << BigInt((x.toString(2).length + 1) >> 1);
  for (;;) {
    const y = (r + x / r) >> 1n;
    if (y >= r) return r;
    r = y;
  }
}

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

/** The volume of the unit n-ball (floats). */
export function ndUnitBallVolume(n: number): number {
  const v = BALL_VOLUME[n];
  return (Number(v.num) / Number(v.den)) * Math.pow(Math.PI, v.k);
}

/** The volume of a ball of radius 1/2 in n dimensions. */
export const ndBallVolume = (n: number): number => ndUnitBallVolume(n) / 2 ** n;

/** The n-volume of the boundary of the container at unit edge (floats), n = d - 1. */
export function ndBoundaryVolume(family: NDFamily, d: number): number {
  const n = d - 1;
  const nf = Number(factorial(n));
  const simplex = Math.sqrt(n + 1) / (nf * 2 ** (n / 2));
  if (family === "simplex") return (d + 1) * simplex;
  if (family === "orthoplex") return 2 ** d * simplex;
  const demi = (1 - 2 ** (n - 1) / nf) / 2 ** (n / 2);
  return 2 * d * demi + 2 ** (d - 1) * simplex;
}

/** The volume bound as a rational Q with s^(2n) >= Q (so s >= Q^(1/2n)), pi and sqrt d rounded the safe way. */
function volumeQ(family: NDFamily, d: number, P: number): { num: bigint; den: bigint } {
  const n = d - 1;
  const v = BALL_VOLUME[n];
  const nf = factorial(n);
  // B = P (num/den) pi^k n!  with pi rounded down:  Bn / Bd
  const Bn = BigInt(P) * v.num * PI_60 ** BigInt(v.k) * nf;
  const Bd = v.den * TEN60 ** BigInt(v.k);
  if (family === "simplex" || family === "orthoplex") {
    // s^n >= B / (2^(n/2) (d+1) sqrt d)         (simplex)
    // s^n >= B / (2^(n/2) 2^d sqrt d)           (orthoplex)
    // squared:  s^(2n) >= B^2 / (2^n c^2 d),  c = d + 1 or 2^d
    const c = family === "simplex" ? BigInt(d + 1) : 2n ** BigInt(d);
    return { num: Bn * Bn, den: Bd * Bd * 2n ** BigInt(n) * c * c * BigInt(d) };
  }
  // demicube: s^n >= B / (2^(n/2) (alpha + beta sqrt d)),  alpha = 2d (n! - 2^(n-1)), beta = 2^(d-1); sqrt d rounded up
  const alpha = BigInt(2 * d) * (nf - 2n ** BigInt(n - 1));
  const beta = 2n ** BigInt(d - 1);
  const sqrtUp = isqrt(BigInt(d) * TEN60 * TEN60) + 1n; // sqrt(d) <= sqrtUp / 10^60
  const denom = alpha * TEN60 + beta * sqrtUp; // (alpha + beta sqrt d) 10^60, rounded up
  return { num: Bn * Bn * TEN60 * TEN60, den: Bd * Bd * 2n ** BigInt(n) * denom * denom };
}

/** The one-ball bound as a rational Q with s^2 >= Q. */
function oneBallQ(family: NDFamily, d: number): { num: bigint; den: bigint } {
  if (family === "simplex" || family === "orthoplex") return { num: BigInt(d * (d - 2)), den: 4n };
  const n = d - 1;
  return { num: BigInt(2 * (n - 1)), den: BigInt(n * n - 5 * n + 8) };
}

/** The exact lower bound on s of P balls on the boundary of the container, as the largest binary64 not above it. */
export function ndBound(family: NDFamily, d: number, P: number): number {
  const n = d - 1;
  const one = oneBallQ(family, d);
  const oneB = largestFloat((x) => compare(x, 2, one.num, one.den) <= 0, Math.sqrt(Number(one.num) / Number(one.den)));
  const q = volumeQ(family, d, P);
  const approx = Math.pow((P * ndBallVolume(n)) / ndBoundaryVolume(family, d), 1 / n);
  const volB = largestFloat((x) => compare(x, 2 * n, q.num, q.den) <= 0, approx);
  return Math.max(oneB, volB);
}

/** Which of the two bounds is the larger (for the listing's label). */
export function ndBoundKind(family: NDFamily, d: number, P: number): "volume" | "one ball" {
  const n = d - 1;
  const one = oneBallQ(family, d);
  const oneB = Math.sqrt(Number(one.num) / Number(one.den));
  const volB = Math.pow((P * ndBallVolume(n)) / ndBoundaryVolume(family, d), 1 / n);
  return volB > oneB ? "volume" : "one ball";
}

/** The bound of a problem by its container name. */
export function ndProblemBound(polytope: string, P: number): number {
  const p = ndParse(polytope);
  if (!p) throw new Error(`ndProblemBound: ${polytope} is none of the simplex, orthoplex or demicube containers`);
  return ndBound(p.family, p.d, P);
}

/** The cell counts and sizes the docs table lists: facets, ridges and hinges of each container. */
export function ndCounts(family: NDFamily, d: number): { facets: number; ridges: number; simplexFacets: number; otherFacets: number } {
  if (family === "simplex") return { facets: d + 1, ridges: ((d + 1) * d) / 2, simplexFacets: d + 1, otherFacets: 0 };
  if (family === "orthoplex") return { facets: 2 ** d, ridges: 2 ** d * (d / 2), simplexFacets: 2 ** d, otherFacets: 0 };
  // demicube: 2 d demicube facets, 2^(d-1) simplex facets; each simplex facet has d ridges (all shared with demicube facets),
  // each demicube facet 2 (d-1) ridges shared with demicube facets and 2^(d-2) with simplex facets
  const dem = 2 * d;
  const sim = 2 ** (d - 1);
  return { facets: dem + sim, ridges: sim * d + (dem * 2 * (d - 1)) / 2, simplexFacets: sim, otherFacets: dem };
}

/** The edge from which a ball at a cell's centre is valid (for the cell kinds that exist), and for pairs of adjacent centres. */
export function ndCentreEdge(polytope: NDName): number {
  const p = ndParse(polytope)!;
  if (p.family === "demicube") return p.d === 5 ? Math.sqrt(6) / 2 : 1;
  return Math.sqrt(p.d * (p.d - 2)) / 2;
}
export function ndPairEdge(polytope: NDName): number {
  const p = ndParse(polytope)!;
  if (p.family === "demicube") return Math.SQRT2;
  return Math.sqrt((p.d * (p.d - 1)) / 2);
}

/**
 * Estimated milliseconds of the server's check per ball ([CALC] BOUNDARY.md 4C.7, 2026-11-02, one core): set from measurements.
 */
export const ND_BALL_MS: Record<NDName, number> = Object.fromEntries(ND_NAMES.map((n) => [n, 5])) as Record<NDName, number>;

export function ndCheckMs(polytope: NDName, P: number): number {
  return ND_BALL_MS[polytope] * P;
}
