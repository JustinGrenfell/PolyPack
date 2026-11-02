/**
 * Structured starts of the d-cube kernel (kernel.ts, capability "seeds"): the cell layouts and the ring cuts.
 *
 * Cell layout. Each of the 2d cells holds floor(P / 2d) or ceil(P / 2d) = m pieces (P mod 2d random cells hold one more), at
 * points of one local layout:
 *  - balls (diameter 1): one at the cell's centre (s = 1); two on a random main diagonal of the inner cube of edge s - 1,
 *    1 apart (s = 1 + 1 / sqrt n); more on the k^n grid of spacing 1, k = ceil(m^(1/n)) (s = k);
 *  - cubes: the k^n grid of unit cubes (s = k, with 1e-6 of room).
 * Every piece lies in the inner cube of its cell, at distance >= 1/2 from every face, so any two pieces are >= 1 apart
 * along any path that crosses a face (BOUNDARY.md 4B): the layout is valid at its edge.
 *
 * Ring cut (a cut of the lattice Z^n by a flat strip). The cells (+-e_a) and (+-e_b) form a ring: the boundary of the square
 * [-h, h]^2 in the axes a, b times [-h, h]^(d-2) is a flat strip of length 4 s (a cyclic row of four cells with trivial
 * holonomy). Pieces sit at m equally spaced points of the ring from a random phase, times a g^(d-2) grid across it (spacing
 * 1, inside [-(h-1/2), h-1/2]); the ring edge is the least s >= g with 4 s >= m, so that neighbours along it are >= 1
 * apart. Some pieces straddle a ridge of the strip. A ball layout is kept only when the validator accepts it.
 */
import type { PieceName } from "../types";
import type { State } from "../search/model";
import { randInt } from "../search/random";
import { axesOf } from "./lattice";

/** Room the cube layouts get over their exact edge. */
const ROOM = 1 + 1e-6;

/** How many pieces each of the 2d cells holds. */
function perCell(d: number, P: number, rand: () => number): number[] {
  const C = 2 * d;
  const counts = new Array<number>(C).fill(Math.floor(P / C));
  const cells = Array.from({ length: C }, (_, i) => i);
  for (let i = 0; i < P % C; i++) {
    const j = i + randInt(rand, C - i);
    [cells[i], cells[j]] = [cells[j], cells[i]];
    counts[cells[i]]++;
  }
  return counts;
}

/** `count` distinct points of the k^n grid of spacing `step` centred at the origin, in random order. */
function gridPoints(n: number, k: number, step: number, count: number, rand: () => number): number[][] {
  const total = k ** n;
  const picked = new Set<number>();
  while (picked.size < count) picked.add(randInt(rand, total));
  return [...picked].map((idx) => {
    const p: number[] = [];
    let r = idx;
    for (let a = 0; a < n; a++) {
      p.push(((r % k) - (k - 1) / 2) * step);
      r = Math.floor(r / k);
    }
    return p;
  });
}

/** The cell layout of P pieces on the d-cube (header). */
export function cellLayout(d: number, piece: PieceName, P: number, rand: () => number): { st: State; s: number } {
  const n = d - 1;
  const C = 2 * d;
  const m = Math.ceil(P / C);
  const k = Math.max(1, Math.ceil(Math.pow(m, 1 / n) - 1e-12));
  let s: number;
  let local: (count: number) => number[][];
  if (piece === "cube") {
    s = k * ROOM;
    local = (count) => gridPoints(n, k, s / k, count, rand);
  } else if (m === 1) {
    s = 1;
    local = () => [new Array<number>(n).fill(0)];
  } else if (m === 2) {
    s = 1 + 1 / Math.sqrt(n);
    local = () => {
      const a = (s - 1) / 2;
      const p = Array.from({ length: n }, () => (rand() < 0.5 ? -a : a));
      return [p, p.map((x) => -x)];
    };
  } else {
    s = k;
    local = (count) => gridPoints(n, k, 1, count, rand);
  }
  const counts = perCell(d, P, rand);
  const cells = new Int32Array(P);
  const U = new Float64Array(n * P);
  let q = 0;
  for (let c = 0; c < C; c++) {
    if (!counts[c]) continue;
    const pts = local(counts[c]);
    for (let i = 0; i < counts[c]; i++) {
      cells[q] = c;
      U.set(pts[i], n * q);
      q++;
    }
  }
  return { st: { cells, U }, s };
}

/** The ring cut (header) of P pieces, or null when the layout would not fit with a cross grid of at most 3 per axis. */
export function ringLayout(
  d: number,
  piece: PieceName,
  P: number,
  rand: () => number,
): { st: State; s: number } | null {
  const n = d - 1;
  // two random axes a < b of R^d
  const a = randInt(rand, d);
  let b = randInt(rand, d - 1);
  if (b >= a) b++;
  const [ra, rb] = a < b ? [a, b] : [b, a];
  const rest = Array.from({ length: d }, (_, i) => i).filter((i) => i !== ra && i !== rb);
  let best: { g: number; s: number } | null = null;
  for (let g = 1; g <= 3; g++) {
    const layers = g ** rest.length;
    const m = Math.ceil(P / layers);
    const s = Math.max(g, m / 4);
    if (!best || s < best.s) best = { g, s };
  }
  if (!best) return null;
  const g = best.g;
  const s = best.s * (piece === "cube" ? ROOM : 1 + 1e-9);
  const h = s / 2;
  const layers = g ** rest.length;
  const m = Math.ceil(P / layers);
  const phase = rand() * 4 * s;
  const cells = new Int32Array(P);
  const U = new Float64Array(n * P);
  const cross = gridPoints(rest.length, g, g > 1 ? (s - 1) / (g - 1) : 1, Math.min(layers, P), () => rand());
  let q = 0;
  for (let j = 0; j < m && q < P; j++) {
    const t = (phase + (j * 4 * s) / m) % (4 * s);
    const tt = (t + h) % (4 * s);
    const seg = Math.floor(tt / s);
    const off = tt - seg * s - h;
    // the R^d point: the ring coordinates and the cross coordinates
    for (const cp of cross) {
      if (q >= P) break;
      const x = new Float64Array(d);
      rest.forEach((ax, i) => (x[ax] = cp[i]));
      // seg 0: cell (ra, +), moving +rb; seg 1: (rb, +), moving -ra; seg 2: (ra, -), moving -rb; seg 3: (rb, -), moving +ra
      let c: number;
      if (seg === 0) {
        x[ra] = h;
        x[rb] = off;
        c = 2 * ra + 1;
      } else if (seg === 1) {
        x[rb] = h;
        x[ra] = -off;
        c = 2 * rb + 1;
      } else if (seg === 2) {
        x[ra] = -h;
        x[rb] = -off;
        c = 2 * ra;
      } else {
        x[rb] = -h;
        x[ra] = off;
        c = 2 * rb;
      }
      const ax = axesOf(d, c >> 1);
      cells[q] = c;
      for (let k = 0; k < n; k++) U[n * q + k] = x[ax[k]];
      q++;
    }
  }
  return q === P ? { st: { cells, U }, s } : null;
}
