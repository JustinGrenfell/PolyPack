/**
 * Moving a centre on the boundary: the exponential map, in ambient coordinates (the nD twin of ../chains.ts traceGeodesic and
 * search/chart.ts for the 4-polytopes).
 *
 * A point x of cell c (its hyperplane H_c) moved along the straight direction d (tangent to H_c) runs inside the cell until it reaches
 * a facet plane n_k . x = s h_k; there it crosses the ridge into the neighbour g = nbr[k]: the point stays, and the direction is carried
 * by the unfolding reflection of the ridge (chains.ts: m = n_c + n_g), an involution that maps H_c onto H_g, so the straight line
 * beyond the ridge in H_c becomes the straight line in H_g. A segment that reaches a hinge (two facet planes at once) has no
 * continuation: the boundary is not flat there. Every crossing is a reflection of the direction, so the linear part of the whole map
 * is a product of reflections, returned as a D x D matrix M (d_end = M d_start) when asked.
 */
import type { NDContainer } from "./container";

/** Crossings after which a walk gives up. */
export const ND_MAX_CROSSINGS = 4096;

export interface NDWalk {
  /** The end cell and point (ambient, edge s), the cells passed, and whether a hinge stopped the walk (the point is then where it did). */
  cell: number;
  x: Float64Array;
  hit: boolean;
  cells: number[];
}

/**
 * Walks from x (cell c, edge s) along d (a vector tangent to H_c, any length) for the length |d|. `M`, when given (D*D numbers),
 * receives the linear part of the map d_start -> d_end (identity when nothing was crossed).
 */
export function walkND(
  C: NDContainer,
  s: number,
  c: number,
  x0: ArrayLike<number>,
  d0: ArrayLike<number>,
  M: Float64Array | null = null,
  hingeTol = 1e-12,
): NDWalk {
  const D = C.D;
  const x = Float64Array.from(x0);
  const d = Float64Array.from(d0);
  let len = 0;
  for (let j = 0; j < D; j++) len += d[j] * d[j];
  len = Math.sqrt(len);
  if (M) {
    M.fill(0);
    for (let j = 0; j < D; j++) M[j * D + j] = 1;
  }
  if (!(len > 0)) return { cell: c, x, hit: false, cells: [c] };
  for (let j = 0; j < D; j++) d[j] /= len;
  let rem = len;
  let cur = c;
  const cells = [c];
  const m = new Float64Array(D);
  for (let step = 0; step < ND_MAX_CROSSINGS; step++) {
    const cl = C.cell(cur);
    const F = cl.nbr.length;
    let best = Infinity;
    let kb = -1;
    for (let k = 0; k < F; k++) {
      let nd = 0;
      let nx = 0;
      for (let j = 0; j < D; j++) {
        nd += cl.nn[k * D + j] * d[j];
        nx += cl.nn[k * D + j] * x[j];
      }
      if (nd <= 1e-15) continue;
      const t = (s * cl.hn[k] - nx) / nd;
      if (t < best) {
        best = t;
        kb = k;
      }
    }
    if (best >= rem || kb < 0) {
      for (let j = 0; j < D; j++) x[j] += rem * d[j];
      return { cell: cur, x, hit: false, cells };
    }
    best = Math.max(best, 0);
    for (let j = 0; j < D; j++) x[j] += best * d[j];
    for (let k = 0; k < F; k++) {
      if (k === kb) continue;
      let nx = 0;
      for (let j = 0; j < D; j++) nx += cl.nn[k * D + j] * x[j];
      if (s * cl.hn[k] - nx < hingeTol * Math.max(1, s)) return { cell: cur, x, hit: true, cells };
    }
    let mm = 0;
    for (let j = 0; j < D; j++) {
      m[j] = cl.n[j] + cl.nn[kb * D + j];
      mm += m[j] * m[j];
    }
    const f = (2 * dotv(m, d, D)) / mm;
    for (let j = 0; j < D; j++) d[j] -= f * m[j];
    if (M) {
      // M <- R M : each column reflected
      for (let col = 0; col < D; col++) {
        let md = 0;
        for (let r = 0; r < D; r++) md += m[r] * M[r * D + col];
        const g = (2 * md) / mm;
        for (let r = 0; r < D; r++) M[r * D + col] -= g * m[r];
      }
    }
    cur = cl.nbr[kb];
    cells.push(cur);
    rem -= best;
  }
  throw new Error("walkND: too many crossings");
}

function dotv(a: Float64Array, b: Float64Array, D: number): number {
  let s = 0;
  for (let j = 0; j < D; j++) s += a[j] * b[j];
  return s;
}
