/**
 * Polytope data from formulas: the simplex, hypercube, orthoplex and demicube of every dimension, and the 4_21 polytope of
 * 8D, computed when first needed instead of shipped as tables (the tables of these were 1.9 MB of source, and a bundled
 * client carries every module it can import). The values agree with the generated exact data (python/ndim/gen_catalogue.py,
 * exact numbers rounded to binary64) to a few ulps: packages/ndim/test/generate.test.ts compares them with that data.
 * Vertex order is the generator's; the order of the facets is not (nothing depends on it).
 *
 * Unit edge, centred at the vertex centroid:
 *   simplex   vertices from the orthonormalised vectors (1^k, -k, 0, ...), scaled by 1/sqrt(2); facet j opposite vertex j.
 *   cube      (+-1/2)^d; facets +-e_i at distance 1/2.
 *   orthoplex +-e_i / sqrt(2); facets (+-1)^d / sqrt(d) at distance 1 / sqrt(2d).
 *   demicube  (+-a)^d with an even number of minus signs, a = 1 / (2 sqrt(2)); simplex facets s / sqrt(d) (an odd number of
 *             minus signs in s) at distance a (d - 2) / sqrt(d), and the cube facets +-e_i at distance a.
 *   4_21      the 240 roots of E8 over sqrt(2); its 19,440 facets are the orbits of two normals under the reflections of E8
 *             (17,280 simplex facets at distance 3/4 and 2,160 orthoplex facets at distance 1/sqrt(2)).
 */
import type { PolytopeData } from "./polytope";

type Data = Pick<PolytopeData, "V" | "N" | "b">;

function simplex(d: number): Data {
  const V: number[] = [];
  const N: number[] = [];
  const R = Math.sqrt(d / (2 * (d + 1)));
  for (let j = 0; j <= d; j++) {
    for (let k = 1; k <= d; k++) {
      const bk = j < k ? 1 : j === k ? -k : 0;
      V.push(bk / Math.sqrt(2 * k * (k + 1)));
    }
  }
  for (let j = 0; j <= d; j++) for (let k = 0; k < d; k++) N.push(-V[j * d + k] / R);
  return { V, N, b: new Array(d + 1).fill(R / d) };
}

/** Sign vectors of length d in binary order with + before -, optionally only those with an even or odd number of minus signs. */
function signVectors(d: number, parity?: 0 | 1): number[][] {
  const out: number[][] = [];
  for (let m = 0; m < 1 << d; m++) {
    const s: number[] = [];
    let minus = 0;
    for (let i = 0; i < d; i++) {
      const neg = (m >> (d - 1 - i)) & 1;
      minus += neg;
      s.push(neg ? -1 : 1);
    }
    if (parity === undefined || minus % 2 === parity) out.push(s);
  }
  return out;
}

const axisFacets = (d: number, dist: number, N: number[], b: number[]) => {
  for (let i = 0; i < d; i++)
    for (const sg of [1, -1]) {
      for (let k = 0; k < d; k++) N.push(k === i ? sg : 0);
      b.push(dist);
    }
};

function cube(d: number): Data {
  const V = signVectors(d).flatMap((s) => s.map((x) => x / 2));
  const N: number[] = [];
  const b: number[] = [];
  axisFacets(d, 0.5, N, b);
  return { V, N, b };
}

function ortho(d: number): Data {
  const V: number[] = [];
  for (let i = 0; i < d; i++)
    for (const sg of [1, -1]) for (let k = 0; k < d; k++) V.push(k === i ? sg * Math.SQRT1_2 : 0);
  const N = signVectors(d).flatMap((s) => s.map((x) => x / Math.sqrt(d)));
  return { V, N, b: new Array(1 << d).fill(1 / Math.sqrt(2 * d)) };
}

function demi(d: number): Data {
  const a = 1 / (2 * Math.SQRT2);
  const V = signVectors(d, 0).flatMap((s) => s.map((x) => x * a));
  const N = signVectors(d, 1).flatMap((s) => s.map((x) => x / Math.sqrt(d)));
  const b = new Array(N.length / d).fill((a * (d - 2)) / Math.sqrt(d));
  axisFacets(d, a, N, b);
  return { V, N, b };
}

/** The 240 roots of E8, scaled to unit vectors: (+-1, +-1, 0^6) and (+-1/2)^8 with an even number of minus signs. */
function e8Roots(): number[][] {
  const R: number[][] = [];
  for (let i = 0; i < 8; i++)
    for (let j = i + 1; j < 8; j++)
      for (const si of [1, -1])
        for (const sj of [1, -1]) {
          const v = new Array(8).fill(0);
          v[i] = si;
          v[j] = sj;
          R.push(v);
        }
  for (const s of signVectors(8, 0)) R.push(s.map((x) => x / 2));
  return R;
}

/**
 * The orthonormal frame the generator puts E8 in: Gram-Schmidt of the roots in the order of e8Roots (the first root becomes
 * the first axis, and so on), so that vertex k of the polytope is the root's coordinates along these axes over sqrt(2).
 */
function e8Frame(roots: number[][]): number[][] {
  const frame: number[][] = [];
  for (const r of roots) {
    const w = r.slice();
    for (const u of frame) {
      const f = w.reduce((a, x, k) => a + x * u[k], 0);
      for (let k = 0; k < 8; k++) w[k] -= f * u[k];
    }
    const len = Math.hypot(...w);
    if (len > 1e-6) frame.push(w.map((x) => x / len));
    if (frame.length === 8) break;
  }
  return frame;
}

function e8Cell(): Data {
  const roots = e8Roots();
  const frame = e8Frame(roots);
  const inFrame = (x: number[]) => frame.map((u) => x.reduce((a, c, k) => a + c * u[k], 0));
  const V = roots.flatMap((r) => inFrame(r.map((x) => x * Math.SQRT1_2)));
  // simple roots (Bourbaki), unit length: reflecting in them generates the Weyl group, which preserves the polytope
  const simple = [
    [0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, 0.5],
    [1, 1, 0, 0, 0, 0, 0, 0],
    [-1, 1, 0, 0, 0, 0, 0, 0],
    [0, -1, 1, 0, 0, 0, 0, 0],
    [0, 0, -1, 1, 0, 0, 0, 0],
    [0, 0, 0, -1, 1, 0, 0, 0],
    [0, 0, 0, 0, -1, 1, 0, 0],
    [0, 0, 0, 0, 0, -1, 1, 0],
  ].map((r) => inFrame(r.map((x) => x * Math.SQRT1_2)));
  const s = Math.SQRT2 / 8;
  const reps: { n: number[]; b: number; count: number }[] = [
    { n: [-0.75, -0.5, -s, -s, -s, -s, -s, -s], b: 0.75, count: 17280 },
    { n: [-Math.SQRT1_2, -Math.SQRT1_2 / 2, -0.25, -0.25, -0.25, -0.25, -0.25, 0.25], b: Math.SQRT1_2, count: 2160 },
  ];
  const N: number[] = [];
  const b: number[] = [];
  for (const rep of reps) {
    const key = (x: number[]) => x.map((c) => Math.round(c * 1e7)).join(",");
    const seen = new Set<string>([key(rep.n)]);
    const queue: number[][] = [rep.n];
    for (let q = 0; q < queue.length; q++) {
      const x = queue[q];
      for (const u of simple) {
        let dot = 0;
        for (let k = 0; k < 8; k++) dot += x[k] * u[k];
        const y = x.map((c, k) => c - 2 * dot * u[k]);
        const len = Math.hypot(...y);
        const z = y.map((c) => c / len);
        const kz = key(z);
        if (!seen.has(kz)) {
          seen.add(kz);
          queue.push(z);
        }
      }
    }
    if (queue.length !== rep.count) throw new Error(`4_21: ${queue.length} facets in an orbit, expected ${rep.count}`);
    for (const x of queue) {
      N.push(...x);
      b.push(rep.b);
    }
  }
  return { V, N, b };
}

/** The formula families (simplex, cube, ortho, demi) of any dimension, and 8:421. */
export const FORMULA_NAMES: readonly string[] = ["simplex", "cube", "ortho", "demi", "421"];

/** The data of a polytope computed from formulas, or null when the code is not one of those (see the file comment). */
export function generatePolytope(code: string): Data | null {
  const m = /^(\d+):(simplex|cube|ortho|demi|421)$/.exec(code);
  if (!m) return null;
  const d = +m[1];
  switch (m[2]) {
    case "simplex":
      return simplex(d);
    case "cube":
      return cube(d);
    case "ortho":
      return ortho(d);
    case "demi":
      return d >= 5 ? demi(d) : null;
    default:
      return d === 8 ? e8Cell() : null;
  }
}
