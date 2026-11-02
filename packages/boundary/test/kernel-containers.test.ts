/**
 * The boundary kernel on the five containers beyond the tesseract (balls only): the generic methods' contract, the exact gradient
 * of the chain metric (in every displacement and the size, across ridges), moves along an independent R^4 trace, random poses
 * along geodesics through the chains, structured seeds, finish and the server's rule, and the test packings.
 */
import { readFileSync } from "node:fs";
import { fileURLToPath } from "node:url";
import { type KernelState, mulberry32, squeeze } from "@polypack/methods";
import { checkKernel } from "@polypack/methods/contract";
import { describe, expect, it } from "vitest";
import type { BoundaryConfig, PolytopeName } from "../src/index";
import {
  ballBounds,
  boundaryKernel,
  geodesicDistance,
  getPolytope,
  traceGeodesic,
  validateBoundary,
} from "../src/index";
import { displacementTo } from "../src/kernel-pose";
import { randomCellPoint } from "../src/search/chart";
import { trace4d } from "./helpers";

const NAMES: PolytopeName[] = ["5cell", "16cell", "24cell", "600cell", "120cell"];
const problem = (polytope: PolytopeName, P: number) => ({ polytope, piece: "ball" as const, P });
const kernel = (polytope: PolytopeName, P: number) => boundaryKernel(problem(polytope, P));

interface Packing {
  id: string;
  polytope: PolytopeName;
  P: number;
  score_cert: number;
  config: BoundaryConfig;
}
const packings = (): Packing[] =>
  JSON.parse(
    readFileSync(
      fileURLToPath(new URL("../../../python/boundary/fixtures/general_packings.json", import.meta.url)),
      "utf8",
    ),
  ).packings;

/** A random start shrunk below its edge with every piece pushed up to 1.5 s, so that overlaps abound and paths cross ridges. */
function displaced(k: ReturnType<typeof kernel>, rng: () => number): KernelState {
  const s = k.randomStart(rng);
  k.setSize(s, s.size * (0.7 + 0.3 * rng()));
  for (let i = 0; i < k.P; i++) k.jiggle(s, i, (1.5 * s.size) / k.steps[0], rng);
  return s;
}

describe("boundary kernel on the other containers: the contract", () => {
  for (const name of NAMES)
    it(`${name} passes checkKernel`, () => {
      const Ps = name === "120cell" || name === "600cell" ? [1, 3] : [1, 4, 9];
      for (const P of Ps) expect(checkKernel(() => kernel(name, P), { seed: P, states: 2 })).toEqual([]);
    });

  it("has no kernel for a cube or a P above the server's largest", () => {
    expect(() => boundaryKernel({ polytope: "24cell", piece: "cube", P: 3 } as never)).not.toThrow;
    expect(() => kernel("24cell", 65)).toThrow(/above the server's largest, 64/);
  });
});

describe("boundary kernel on the other containers: energy", () => {
  it("has the gradient of its weighted energy, in every displacement and the size", () => {
    for (const name of NAMES) {
      const P = name === "120cell" || name === "600cell" ? 2 : 3;
      const k = kernel(name, P);
      const rng = mulberry32(P + 11);
      const n = 3 * P;
      let tested = 0;
      for (let t = 0; t < 8 && tested < 3; t++) {
        const s = displaced(k, rng);
        const w = Float64Array.from({ length: P }, () => 0.5 + rng());
        const g = new Float64Array(n + 1);
        if (!(k.energy(s, g, w) > 1e-8)) continue;
        tested++;
        for (let i = 0; i <= n; i++) {
          const h = 1e-7 * (i === n ? s.size : 1);
          const p = k.clone(s);
          const m = k.clone(s);
          if (i === n) {
            p.size += h;
            m.size -= h;
          } else {
            p.x[i] += h;
            m.x[i] -= h;
          }
          const fd = (k.energy(p, undefined, w) - k.energy(m, undefined, w)) / (2 * h);
          expect(Math.abs(fd - g[i]), `${name} ${i}`).toBeLessThan(2e-5 * Math.max(1, Math.abs(fd)));
        }
      }
      expect(tested, name).toBeGreaterThan(0);
    }
  });

  it("is linear in the weights; a piece's energy is its own terms and its pairs", () => {
    for (const name of ["24cell", "5cell"] as const) {
      const P = 4;
      const k = kernel(name, P);
      const rng = mulberry32(31);
      for (let t = 0; t < 3; t++) {
        const s = displaced(k, rng);
        const E = k.energy(s);
        let sum = 0;
        for (let i = 0; i < P; i++)
          sum += k.energy(
            s,
            undefined,
            Float64Array.from({ length: P }, (_, j) => (j === i ? 1 : 0)),
          );
        expect(Math.abs(sum - E)).toBeLessThanOrEqual(1e-9 * Math.max(1, E));
        for (let i = 0; i < P; i++) expect(k.clear(s, i, 0)).toBe(k.pieceEnergy(s, i) === 0);
      }
    }
  });
});

describe("boundary kernel on the other containers: charts and poses", () => {
  it("moves a piece along the straight path an independent R^4 trace follows", () => {
    const rng = mulberry32(5);
    for (const name of NAMES) {
      const poly = getPolytope(name);
      const k = kernel(name, 1);
      let checked = 0;
      for (let t = 0; t < 40; t++) {
        const s = 1.5 + rng();
        const cell = Math.floor(rng() * poly.cells.length);
        const u = randomCellPoint(poly, cell, s, rng, 0.5);
        const dir = [0, 1, 2].map(() => 2 * rng() - 1);
        const norm = Math.hypot(...dir);
        const L = s * (0.3 + 3 * rng());
        const ref = trace4d(poly, s, cell, u, dir, L);
        if (!ref) continue;
        const st = k.fromConfig({ s, pieces: [{ cell, u }] })!;
        for (let m = 0; m < 3; m++) st.x[m] = (dir[m] / norm) * L;
        const p = (k.toConfig(st) as BoundaryConfig).pieces[0];
        expect(p.cell).toBe(ref.cell);
        for (let m = 0; m < 3; m++) expect(Math.abs(p.u[m] - ref.u[m])).toBeLessThan(1e-8);
        checked++;
      }
      expect(checked, name).toBeGreaterThan(20);
    }
  });

  it("develops a target through the chains along the shortest geodesic: the walk ends on the target", () => {
    const rng = mulberry32(9);
    for (const name of NAMES) {
      const poly = getPolytope(name);
      let found = 0;
      let missed = 0;
      for (let t = 0; t < 30; t++) {
        // an edge from the one-ball bound up to 2.5 times it, where packings live
        const s = ballBounds(name, 1).value * (1 + 1.5 * rng());
        const c0 = Math.floor(rng() * poly.cells.length);
        const a = randomCellPoint(poly, c0, s, rng, 0.7);
        const c1 = Math.floor(rng() * poly.cells.length);
        const b = randomCellPoint(poly, c1, s, rng, 0.7);
        const v = displacementTo(c0, a, c1, b, s, poly);
        if (!v) {
          missed++;
          continue;
        }
        found++;
        const len = Math.hypot(...v);
        const walk = traceGeodesic(poly, s, c0, a, v, len);
        expect(walk.hitHinge).toBe(false);
        expect(walk.cell).toBe(c1);
        for (let m = 0; m < 3; m++) expect(Math.abs(walk.u[m] - b[m])).toBeLessThan(1e-8);
        // it is the shortest geodesic
        const d = geodesicDistance(poly, s, { cell: c0, u: a }, { cell: c1, u: b }, len * 1.05 + 1e-9);
        expect(len).toBeLessThanOrEqual(d * (1 + 1e-9));
      }
      // a far target of a big container may need more chains than POSE_CHAINS and is given up (the kernel then walks a geodesic):
      // almost always on the 120-cell and 600-cell, whose cells are many and small beside the distances
      expect(found + missed).toBe(30);
      if (name !== "600cell" && name !== "120cell")
        expect(found, `${name} found ${found}, missed ${missed}`).toBeGreaterThan(10);
    }
  });

  it("finds every target within 1.2 of the anchor, at the walk's own length or less", () => {
    const rng = mulberry32(19);
    for (const name of NAMES) {
      const poly = getPolytope(name);
      let tried = 0;
      for (let t = 0; t < 40 && tried < 15; t++) {
        const s = 1.45 + 1.5 * rng();
        const c0 = Math.floor(rng() * poly.cells.length);
        const a = randomCellPoint(poly, c0, s, rng, 0.7);
        const dir = [0, 1, 2].map(() => 2 * rng() - 1);
        const L = 0.2 + 1.0 * rng();
        const end = traceGeodesic(poly, s, c0, a, dir, L);
        if (end.hitHinge) continue;
        tried++;
        const v = displacementTo(c0, a, end.cell, end.u, s, poly);
        expect(v, `${name} ${t}`).not.toBeNull();
        expect(Math.hypot(...v!)).toBeLessThanOrEqual(L * (1 + 1e-9));
        const back = traceGeodesic(poly, s, c0, a, v!, Math.hypot(...v!));
        expect(back.cell).toBe(end.cell);
        for (let m = 0; m < 3; m++) expect(Math.abs(back.u[m] - end.u[m])).toBeLessThan(1e-8);
      }
      expect(tried, name).toBeGreaterThan(10);
    }
  });

  it("randomPose changes one piece's displacement and leaves the anchors", () => {
    const rng = mulberry32(2);
    for (const name of NAMES) {
      const k = kernel(name, 3);
      const s = k.randomStart(rng);
      const before = k.clone(s);
      k.randomPose(s, 1, rng);
      expect(s.aux).toEqual(before.aux);
      expect([...s.x.slice(0, 3)]).toEqual([...before.x.slice(0, 3)]);
      expect([...s.x.slice(6, 9)]).toEqual([...before.x.slice(6, 9)]);
      // the piece moves, whichever way the pose was drawn (a target through the chains, or the end of a random geodesic)
      expect([...s.x.slice(3, 6)]).not.toEqual([...before.x.slice(3, 6)]);
    }
  });
});

describe("boundary kernel on the other containers: starts, seeds and finish", () => {
  it("finishes a squeezed random start into a valid configuration above the bound", () => {
    for (const name of NAMES)
      for (const P of name === "120cell" || name === "600cell" ? [2] : [2, 6]) {
        const k = kernel(name, P);
        const s = k.randomStart(mulberry32(100 + P));
        squeeze(k, s);
        const f = k.finish(s)!;
        expect(f, `${name} ${P}`).not.toBeNull();
        expect(f.score).toBeGreaterThanOrEqual(ballBounds(name, P).value * (1 - 1e-9));
        expect(validateBoundary(problem(name, P), f.config, { tol: 0 }).reasons).toEqual([]);
        expect(k.validate(f.config).valid).toBe(true);
      }
  });

  it("seeds valid layouts: balls at the centres of P cells (P up to the cells), or a sequential random layout", () => {
    const rng = mulberry32(4);
    for (const name of NAMES) {
      const C = getPolytope(name).cells.length;
      for (const P of C >= 24 && name !== "120cell"
        ? [3, 9]
        : name === "5cell"
          ? [3, 5, 9]
          : name === "16cell"
            ? [4, 16]
            : [3]) {
        const k = kernel(name, P);
        const seeds = k.seeds!(rng, 2);
        expect(seeds.length).toBe(2);
        for (const q of seeds) {
          expect(k.energy(q)).toBeLessThan(1e-20);
          expect(k.validate(k.toConfig(q)).reasons, `${name} ${P}`).toEqual([]);
        }
      }
    }
    // one ball per cell of the 24-cell: valid from s = sqrt 6 / 2 (adjacent cell centres are 2 r s apart across a ridge)
    expect(kernel("24cell", 24).seeds!(rng, 1)[0].size).toBeCloseTo(Math.sqrt(6) / 2, 6);
  });

  it("reads, converts, validates and finishes the test packings, at or below their score", () => {
    for (const r of packings().filter((p) => p.P <= 6)) {
      const k = kernel(r.polytope, r.P);
      expect(k.validate(r.config)).toEqual({ valid: true, score: r.score_cert, reasons: [] });
      const st = k.fromConfig(r.config)!;
      expect(st.size).toBe(r.config.s);
      expect(k.energy(st)).toBe(0);
      const f = k.finish(st)!;
      expect(f.score).toBeLessThanOrEqual(r.score_cert * (1 + 1e-9));
      expect(validateBoundary(problem(r.polytope, r.P), f.config, { tol: 0 }).reasons).toEqual([]);
    }
  });

  it("validates as the server does: the lower bound and its slack, the cells' range, the server's words", () => {
    const k = kernel("24cell", 1);
    const b = ballBounds("24cell", 1).value;
    const one = { cell: 0, u: [0, 0, 0] };
    expect(k.validate({ s: 1.1, pieces: [one] }).valid).toBe(true);
    expect(k.validate({ s: b * (1 - 1.1e-9), pieces: [one] })).toEqual({
      valid: false,
      score: b * (1 - 1.1e-9),
      reasons: [`the size ${b * (1 - 1.1e-9)} is below the lower bound ${b} of this problem`],
    });
    expect(k.validate({ s: 2, pieces: [{ cell: 24, u: [0, 0, 0] }] })).toEqual({
      valid: false,
      score: null,
      reasons: ["pieces[0].cell must be an integer from 0 to 23"],
    });
    expect(k.validate({ s: 2, pieces: [{ cell: 23, u: [0, 0, 0] }] }).valid).toBe(true);
    expect(kernel("600cell", 1).validate({ s: 2, pieces: [{ cell: 599, u: [0, 0, 0] }] }).valid).toBe(true);
    expect(kernel("5cell", 1).validate({ s: 2, pieces: [{ cell: 5, u: [0, 0, 0] }] }).reasons).toEqual([
      "pieces[0].cell must be an integer from 0 to 4",
    ]);
  });
});
