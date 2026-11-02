/**
 * Structured starts of the boundary kernel (kernel.ts, capability "seeds"): the cell layouts, the
 * nofold packings of docs/spaces/BOUNDARY.md section 6. Each of the 8 cells holds floor(P / 8) or
 * ceil(P / 8) = n pieces (P mod 8 random cells hold one more), at points of one local layout:
 *  - balls (diameter 1): one at the cell's centre (s = 1); two on a random main diagonal of the inner
 *    cube of edge s - 1, 1 apart (s = 1 + 1 / sqrt 3); more on the k^3 grid of spacing 1,
 *    k = ceil(cbrt n) (s = k);
 *  - cubes: the k^3 grid of unit cubes (s = k, with the two rings' hair of room).
 * Every piece lies in the closed cell that anchors it, so the layout is valid at its edge (pieces
 * may touch across a 2-face). The cubes' other structured start, the search's two rings, is the
 * cube model's constructedStart.
 *
 * The other containers (balls only), generalLayout: with P at most the number of cells C, balls at the centres of P
 * random cells; two adjacent centres are 2 r apart across their common facet (r the cell's inradius at edge s), so the
 * layout is valid from s = 1 / (2 r) (and from 1 / (2 rho), the one-ball bound, rho the cell's midradius): every
 * other pair of centres is farther and a loop around an edge is longer. With more balls than cells, a sequential random
 * layout: balls are added at random points, each accepted when it conflicts with none placed (the chain metric's
 * violation of that ball's terms, seen from its own chains), at an edge raised by 6% after every failed round.
 */
import { containerFacts } from "./bounds";
import { ChainMetric } from "./search/balls";
import { randomCellPoint } from "./search/chart";
import type { State } from "./search/model";
import { randInt } from "./search/random";
import type { PieceName, Polytope } from "./types";
import { BALL_RADIUS } from "./units";

/** A ball's diameter, the least distance of two centres. */
const D = 2 * BALL_RADIUS;

/** Room the cube layouts get over their exact edge, as the search's two rings (1e-6 relative). */
const ROOM = 1 + 1e-6;

/** How many pieces each of the 8 cells holds: floor(P / 8), plus one in P mod 8 random cells. */
function perCell(P: number, rand: () => number): number[] {
  const counts = new Array<number>(8).fill(Math.floor(P / 8));
  const cells = [0, 1, 2, 3, 4, 5, 6, 7];
  for (let i = 0; i < P % 8; i++) {
    const j = i + randInt(rand, 8 - i);
    [cells[i], cells[j]] = [cells[j], cells[i]];
    counts[cells[i]]++;
  }
  return counts;
}

/** The nofold layout of P pieces: each cell's pieces at points of one local layout, in random order. */
export function cellLayout(piece: PieceName, P: number, rand: () => number): { st: State; s: number } {
  const n = Math.ceil(P / 8);
  const k = Math.ceil(Math.cbrt(n) - 1e-12);
  let s: number;
  let local: (rand: () => number) => number[][];
  if (piece === "cube") {
    // k^3 unit cubes, spaced s / k >= 1 apart
    s = k * ROOM;
    local = () => grid(k, s / k);
  } else if (n === 1) {
    // one ball per cell, centred: the cell's edges are at s sqrt 2 / 2 >= the radius, so s = 2 r (site units)
    s = D;
    local = () => [[0, 0, 0]];
  } else if (n === 2) {
    // the two centres on a main diagonal of the inner cube of edge s - D, D apart
    s = D + D / Math.sqrt(3);
    local = (r) => {
      const a = s / 2 - BALL_RADIUS;
      const d = [0, 1, 2].map(() => (r() < 0.5 ? -a : a));
      return [d, d.map((x) => -x)];
    };
  } else {
    s = D * k;
    local = () => grid(k, D);
  }
  const counts = perCell(P, rand);
  const cells = new Int32Array(P);
  const U = new Float64Array(3 * P);
  let q = 0;
  for (let c = 0; c < 8; c++) {
    const pts = local(rand);
    // a random choice of counts[c] of the cell's points
    for (let i = 0; i < counts[c]; i++) {
      const j = i + randInt(rand, pts.length - i);
      [pts[i], pts[j]] = [pts[j], pts[i]];
      cells[q] = c;
      U.set(pts[i], 3 * q);
      q++;
    }
  }
  return { st: { cells, U }, s };
}

/** The k^3 points of spacing d centred at the origin. */
function grid(k: number, d: number): number[][] {
  const out: number[][] = [];
  const c = (a: number) => (a - (k - 1) / 2) * d;
  for (let a = 0; a < k; a++) for (let b = 0; b < k; b++) for (let e = 0; e < k; e++) out.push([c(a), c(b), c(e)]);
  return out;
}

/** The edge from which balls at the centres of adjacent cells of this container are valid (header). */
export function centreLayoutEdge(poly: Polytope): number {
  // distance from the centre to a facet plane, and from the centre to the nearest edge, at s = 1, over the cells
  let r = Infinity;
  for (const c of poly.cells) for (const f of c.facets) r = Math.min(r, f.inside.h);
  const rho = containerFacts(poly.name).midradius;
  return Math.max(D / (2 * r), D / (2 * rho)) * (1 + 1e-9);
}

/** Balls on any container: the centre layout when P <= C, else a sequential random layout (header). */
export function generalLayout(poly: Polytope, P: number, rand: () => number): { st: State; s: number } {
  const C = poly.cells.length;
  const cells = new Int32Array(P);
  const U = new Float64Array(3 * P);
  let s = centreLayoutEdge(poly);
  if (P <= C) {
    const order = Array.from({ length: C }, (_, i) => i);
    for (let i = 0; i < P; i++) {
      const j = i + randInt(rand, C - i);
      [order[i], order[j]] = [order[j], order[i]];
      cells[i] = order[i];
    }
    return { st: { cells, U }, s };
  }
  const metric = new ChainMetric(poly);
  // working lengths are twice the site's (units.ts): the metric's balls have radius 1
  const W = 1 / BALL_RADIUS;
  for (let round = 0; round < 60; round++, s *= 1.06) {
    let placed = 0;
    for (let tries = 0; placed < P && tries < 400 * P; tries++) {
      cells[placed] = randInt(rand, C);
      const p = randomCellPoint(poly, cells[placed], s, rand, 0.6);
      for (let k = 0; k < 3; k++) U[3 * placed + k] = p[k];
      const sub = {
        cells: cells.subarray(0, placed + 1),
        U: Float64Array.from(U.subarray(0, 3 * placed + 3), (x) => x * W),
      };
      if (metric.evaluate(sub.cells, sub.U, s * W, null, true, { focus: placed }) <= 0) placed++;
    }
    if (placed === P) return { st: { cells, U }, s };
  }
  // no sequential layout found in 60 rounds: the last, loosest one is returned and the caller's finish() judges it
  return { st: { cells, U }, s };
}
