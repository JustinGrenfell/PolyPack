/**
 * Structured starts of the nD kernel (kernel.ts, capability "seeds"), docs/spaces/BOUNDARY.md 4D.7.
 *
 * Exact layouts (valid at the edge they are built for, derived in 4D.7):
 *  - one ball at a cell's centre: the centre of a simplex or orthoplex cell is at distance s / sqrt(d (d - 2)) from the hinges, so from
 *    s = sqrt(d (d-2)) / 2 (the one-ball bound); the centre of a demicube facet (a (d-1)-demicube) is s/2 from its nearest hinges
 *    for d >= 6, so from s = 1, and from s = sqrt 6 / 2 for d = 5;
 *  - simplex: P <= d + 1 balls at the centres of P cells. Every two cells are adjacent and two adjacent centres are 2 r apart along the straight
 *    path across their ridge, r = s / sqrt(2 n (n + 1)) the inradius of the cell (n = d - 1), so the layout is valid from s = sqrt(d (d-1) / 2);
 *  - orthoplex: P <= 2 balls at the centres of two adjacent cells, from the same s;
 *  - demicube: P <= 2d balls at the centres of 2d demicube facets: two adjacent facets' centres are 2 (s/(2 sqrt 2)) = s / sqrt 2 apart across their
 *    right-angled ridge, so valid from s = sqrt 2 (and the opposite facets are farther).
 * Otherwise a sequential random layout: balls are added one at a time at random points, each accepted when it conflicts with none placed
 * (the model's violation), the edge raised by 6% (the placed balls scaled with it) after every failed round.
 */
import type { State } from "../search/model";
import { randInt } from "../search/random";
import { ndBound, ndCentreEdge, ndPairEdge } from "./bounds";
import { NDBallModel } from "./model";
import { getND } from "./container";
import type { NDName } from "./names";
import { ndParse } from "./names";

/** Room the exact layouts get over their derived edge (relative). */
const ROOM = 1 + 1e-9;

function pick(count: number, from: number, rand: () => number): number[] {
  const pool = Array.from({ length: from }, (_, i) => i);
  for (let i = 0; i < count; i++) {
    const j = i + randInt(rand, from - i);
    [pool[i], pool[j]] = [pool[j], pool[i]];
  }
  return pool.slice(0, count);
}

/** An exact layout of P balls, or null when none is known for this container and P. */
export function ndExactLayout(polytope: NDName, P: number, rand: () => number): { st: State; s: number } | null {
  const C = getND(polytope);
  const n = C.n;
  const mk = (cells: number[], s: number): { st: State; s: number } => ({
    st: { cells: Int32Array.from(cells), U: new Float64Array(n * cells.length) },
    s,
  });
  if (P === 1) return mk([randInt(rand, C.family === "demicube" ? 2 * C.d : C.C)], ndCentreEdge(polytope) * ROOM);
  if (C.family === "simplex" && P <= C.d + 1) return mk(pick(P, C.C, rand), ndPairEdge(polytope) * ROOM);
  if (C.family === "demicube" && P <= 2 * C.d) return mk(pick(P, 2 * C.d, rand), ndPairEdge(polytope) * ROOM);
  if (C.family === "orthoplex" && P === 2) {
    const c = randInt(rand, C.C);
    return mk([c, c ^ (1 << randInt(rand, C.d))], ndPairEdge(polytope) * ROOM);
  }
  return null;
}

/** The sequential random layout (header). */
export function ndSequentialLayout(polytope: NDName, P: number, rand: () => number): { st: State; s: number } {
  const p = ndParse(polytope)!;
  const C = getND(polytope);
  const n = C.n;
  const models = new Map<number, NDBallModel>();
  const model = (k: number) => {
    let m = models.get(k);
    if (!m) models.set(k, (m = new NDBallModel(polytope, k)));
    return m;
  };
  let s = Math.max(ndBound(p.family, p.d, P), ndCentreEdge(polytope)) * 1.1;
  const cells: number[] = [];
  const U: number[][] = [];
  const state = (extra?: { c: number; u: number[] }): State => {
    const k = cells.length + (extra ? 1 : 0);
    const st: State = { cells: new Int32Array(k), U: new Float64Array(n * k) };
    cells.forEach((c, i) => {
      st.cells[i] = c;
      st.U.set(U[i], i * n);
    });
    if (extra) {
      st.cells[k - 1] = extra.c;
      st.U.set(extra.u, (k - 1) * n);
    }
    return st;
  };
  while (cells.length < P) {
    let placed = false;
    for (let t = 0; t < 60 && !placed; t++) {
      const c = model(1).viableCell(rand, s);
      const u = model(1).randomPoint(c, s, rand, 0.4 / Math.sqrt(n));
      const k = cells.length + 1;
      if (model(k).violation(state({ c, u }), s) <= -1e-9 || (k === 1 && model(1).violation(state({ c, u }), s) <= 0)) {
        cells.push(c);
        U.push(u);
        placed = true;
      }
    }
    if (!placed) {
      s *= 1.06;
      for (let i = 0; i < U.length; i++) U[i] = U[i].map((v) => v * 1.06);
    }
  }
  return { st: state(), s };
}
