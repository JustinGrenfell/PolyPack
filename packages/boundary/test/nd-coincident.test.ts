/**
 * Regressions of the C3 review (docs/spaces/BOUNDARY.md 4D.2, 4D.6): two different balls at (nearly) one point overlap, and the ridge snap and
 * the develop prune leave no band in which a neighbour is neither an anchor nor developed into. The exact twin is in python/boundary/test_ndpoly.py.
 */
import { readFileSync } from "node:fs";
import { fileURLToPath } from "node:url";
import { describe, expect, it } from "vitest";
import {
  type NDName,
  ND_MAX_OPS,
  ND_MAX_WORK,
  boundaryKernel,
  getND,
  measureND,
  ndCentreEdge,
  ndMaxOps,
  ndProblemBound,
  ndServerRule,
  validateND,
} from "../src";

type Cfg = { s: number; pieces: { cell: number; u: number[] }[] };
const problem = (polytope: NDName, P: number) => ({ polytope, piece: "ball" as const, P });
const zeros = (n: number) => new Array<number>(n).fill(0);
const gap: { cases: { polytope: NDName; s: number; pieces: Cfg["pieces"]; label: string }[] } = JSON.parse(
  readFileSync(fileURLToPath(new URL("../../../python/boundary/fixtures/nd_gap_cases.json", import.meta.url)), "utf8"),
);

/** Every way the repository judges a configuration: both float tolerances, the server's rule (with its bound) and the strict one. */
function accepted(name: NDName, cfg: Cfg): string[] {
  const P = cfg.pieces.length;
  const out: string[] = [];
  for (const tol of [0, 1e-11, 1e-9]) if (validateND(problem(name, P), cfg, { tol }).ok) out.push(`tol ${tol}`);
  if (ndServerRule(problem(name, P), cfg, 0).valid) out.push("server rule");
  return out;
}

describe("balls at one point", () => {
  const cases: [NDName, number][] = [
    ["5simplex", 2.7],
    ["5simplex", 6],
    ["8simplex", 6],
    ["10simplex", 6],
    ["6orthoplex", 6],
    ["10orthoplex", 4.48],
    ["5demicube", 6],
    ["6demicube", 3],
    ["9demicube", 1.0001],
    ["10demicube", 6],
  ];
  it("refuses 2 and 24 balls stacked at one point, and the bound of the problem does not rescue them", () => {
    for (const [name, s] of cases)
      for (const P of [2, 24]) {
        const cfg: Cfg = { s, pieces: Array.from({ length: P }, () => ({ cell: 0, u: zeros(getND(name).n) })) };
        expect(accepted(name, cfg), `${name} ${P} s=${s}`).toEqual([]);
        const v = validateND(problem(name, P), cfg, { tol: 1e-9 });
        expect(v.reasons[0], name).toMatch(/Balls 0 and 1 overlap/);
        expect(ndServerRule(problem(name, P), cfg, ndProblemBound(name, P)).valid, name).toBe(false);
      }
  });

  it("refuses pairs closer than 1e-10 max(1, s), and reports their distance (offsets 0 ... 2e-10 and 1e-9)", () => {
    for (const name of ["5simplex", "7orthoplex", "6demicube", "10simplex"] as NDName[])
      for (const off of [0, 1e-15, 1e-12, 5e-11, 2e-10, 1e-9, 1e-6]) {
        const s = 6;
        const u2 = zeros(getND(name).n);
        u2[0] = off;
        const cfg: Cfg = {
          s,
          pieces: [
            { cell: 0, u: zeros(u2.length) },
            { cell: 0, u: u2 },
          ],
        };
        expect(accepted(name, cfg), `${name} offset ${off}`).toEqual([]);
        const m = measureND(problem(name, 2), cfg);
        expect(m.pairs).toHaveLength(1);
        expect(m.pairs[0].dist, `${name} offset ${off}`).toBeLessThan(1);
      }
  });

  it("refuses the same point named under two neighbouring cells", () => {
    for (const name of ["5simplex", "6orthoplex", "6demicube", "8simplex"] as NDName[]) {
      const C = getND(name);
      const s = ndCentreEdge(name) * 1.5;
      const c = 0;
      const cl = C.cell(c);
      for (const k of [0, 1]) {
        const g = cl.nbr[k];
        // the ridge's own point (the engine keeps it in the ridge planes), read in both cells
        const rp = Float64Array.from({ length: C.D }, (_, j) => s * cl.rp[k * C.D + j]);
        const cfg: Cfg = {
          s,
          pieces: [
            { cell: c, u: Array.from(C.toLocal(c, rp, s)) },
            { cell: g, u: Array.from(C.toLocal(g, rp, s)) },
          ],
        };
        expect(accepted(name, cfg), `${name} cells ${c}/${g}`).toEqual([]);
      }
    }
  });

  it("still accepts a ball on a ridge alone (its own zero-length loop is not an overlap)", () => {
    const name: NDName = "6orthoplex";
    const C = getND(name);
    const s = ndCentreEdge(name) * 1.5;
    const cl = C.cell(0);
    const rp = Float64Array.from({ length: C.D }, (_, j) => s * cl.rp[0 * C.D + j]);
    const alone: Cfg = {
      s: ndCentreEdge(name) * 2,
      pieces: [{ cell: 0, u: Array.from(C.toLocal(0, rp, ndCentreEdge(name) * 2)) }],
    };
    expect(validateND(problem(name, 1), alone, { tol: 1e-9 }).ok).toBe(true);
    void s;
  });

  it("gives the kernel's energy a penalty for balls at one point", () => {
    for (const name of ["5simplex", "7orthoplex", "6demicube"] as NDName[]) {
      const k = boundaryKernel(problem(name, 2));
      const q = k.fromConfig({
        s: 6,
        pieces: [
          { cell: 0, u: zeros(getND(name).n) },
          { cell: 0, u: zeros(getND(name).n) },
        ],
      } as never);
      expect(k.energy(q!), name).toBeGreaterThan(0.5);
    }
  });
});

describe("the ridge snap and the develop prune leave no gap (4D.6)", () => {
  it("refuses the 22 pairs 0.9 apart across a ridge that fell in the band, at tolerance 0, 1e-11 and 1e-9", () => {
    expect(gap.cases).toHaveLength(22);
    for (const c of gap.cases)
      for (const tol of [0, 1e-11, 1e-9]) {
        const v = validateND(problem(c.polytope, 2), { s: c.s, pieces: c.pieces }, { tol });
        expect(v.ok, `${c.label} ${c.polytope} s=${c.s} tol=${tol}`).toBe(false);
      }
  });
});

describe("the work cap counts operations as well as steps (4D.8)", () => {
  const tie = (name: NDName, P: number, eps: number): Cfg => ({
    s: 1 + eps,
    pieces: Array.from({ length: P }, (_, i) => ({ cell: i, u: zeros(getND(name).n) })),
  });
  it("cuts off a crafted near-tie on the big demicubes by operations, the same on every run, long before 60,000 steps", () => {
    for (const name of ["7demicube", "9demicube", "10demicube"] as NDName[]) {
      const runs = [0, 1].map(() => {
        const stats = { steps: 0, ops: 0 };
        const v = validateND(problem(name, 1), tie(name, 1, 1e-12), { tol: 1e-9, stats });
        return { ok: v.ok, reason: v.reasons[0], stats };
      });
      expect(runs[0].ok, name).toBe(false);
      expect(runs[0].reason, name).toMatch(/work limit \(\d+ operations\)/);
      expect(runs[1]).toEqual(runs[0]);
      expect(runs[0].stats.ops, name).toBeLessThan(ndMaxOps(name) * 1.01);
      expect(runs[0].stats.steps, name).toBeLessThan(ND_MAX_WORK);
    }
  });
  it("accepts the facet-centre balls of the 8-demicube (every hinge tied) within its allowance, and the packings stay far below the cap", () => {
    const stats = { steps: 0, ops: 0 };
    const v = validateND(problem("8demicube", 2), tie("8demicube", 2, 1e-10), { tol: 1e-9, stats });
    expect(v.ok, v.reasons[0]).toBe(true);
    expect(stats.ops).toBeLessThan(ndMaxOps("8demicube"));
    expect(stats.ops).toBeGreaterThan(ND_MAX_OPS);
  });
});
