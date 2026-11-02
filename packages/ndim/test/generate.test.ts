import { describe, expect, it } from "vitest";
import { getPolytope, loadPolytope, polytopeInfo, POLYTOPE_CODES, type PolytopeData } from "../src";
import { EXACT_FAMILIES } from "./fixtures/exact-families";
import { EXACT_421 } from "./fixtures/exact-421";

// The formulas of src/generate.ts against the exact tables of python/ndim/gen_catalogue.py (exact numbers rounded to binary64).
const TOL = 1e-13;

/** Rows of width w, sorted by their values rounded to 1e-8 (the facet order of the formulas is not the generator's). */
function sortedRows(a: ArrayLike<number>, w: number, extra?: ArrayLike<number>): number[][] {
  const rows: number[][] = [];
  for (let i = 0; i < a.length / w; i++) {
    const r = Array.from({ length: w }, (_, k) => a[i * w + k]);
    if (extra) r.push(extra[i]);
    rows.push(r);
  }
  const key = (r: number[]) => r.map((x) => Math.round(x * 1e8) + 0).join(",");
  return rows.sort((p, q) => (key(p) < key(q) ? -1 : key(p) > key(q) ? 1 : 0));
}

function expectSame(got: Float64Array, ref: PolytopeData, d: number) {
  const v = sortedRows(got, d);
  const w = sortedRows(ref.V, d);
  expect(v.length, `${ref.code} vertices`).toBe(w.length);
  let err = 0;
  for (let i = 0; i < v.length; i++) for (let k = 0; k < d; k++) err = Math.max(err, Math.abs(v[i][k] - w[i][k]));
  expect(err, `${ref.code} vertices`).toBeLessThan(TOL);
}

function check(ref: PolytopeData) {
  const info = polytopeInfo(ref.code);
  const P = getPolytope(ref.code);
  expectSame(P.V, ref, info.d);
  expect(P.nF).toBe(ref.b.length);
  const f = sortedRows(P.N, info.d, P.b);
  const g = sortedRows(ref.N, info.d, ref.b);
  expect(f.length, `${ref.code} facets`).toBe(g.length);
  let err = 0;
  for (let i = 0; i < f.length; i++) for (let k = 0; k <= info.d; k++) err = Math.max(err, Math.abs(f[i][k] - g[i][k]));
  expect(err, `${ref.code} facets`).toBeLessThan(TOL);
}

describe("polytopes computed from formulas equal the exact tables", () => {
  it("simplex, hypercube, orthoplex and demicube in 4 to 10 dimensions", () => {
    expect(EXACT_FAMILIES.length).toBe(7 + 7 + 7 + 6);
    for (const ref of EXACT_FAMILIES) check(ref);
  });
  it("4_21 in 8 dimensions", async () => {
    await loadPolytope("8:421");
    check(EXACT_421[0]);
  });
  it("every catalogue entry can be loaded, with the counts of the catalogue", async () => {
    for (const code of POLYTOPE_CODES) {
      const P = await loadPolytope(code);
      const info = polytopeInfo(code);
      expect([P.nV, P.nF]).toEqual([info.nV, info.nF]);
      expect(P.V.length).toBe(info.nV * info.d);
    }
  });
});
