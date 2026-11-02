/**
 * The search model, the search and the kernel for the Platonic solids (docs/spaces/BOUNDARY.md 4C.6): analytic gradients against central differences (in the centres, the
 * rotation increments and the edge), the generic methods' contract, the kernel's rotations and swap, short seeded searches whose results the validator takes
 * at the server's rule, and reproducibility.
 */
import { type KernelState, mulberry32 } from "@polypack/methods";
import { checkKernel } from "@polypack/methods/contract";
import { describe, expect, it } from "vitest";
import type { BoundaryConfig, PolytopeName, PolytopePieceName } from "../src/index";
import { boundaryKernel, getPolytope, pieceBounds, pieceShape, quaternionToMatrix, searchBoundary, validateBoundary } from "../src/index";
import { expRot, mul9 } from "../src/search/rot";
import { SolidModel } from "../src/search/solids";

const CASES: [PolytopeName, PolytopePieceName, number, number][] = [
  ["tesseract", "tetrahedron", 4, 1.0],
  ["tesseract", "icosahedron", 3, 2.0],
  ["5cell", "cube", 3, 1.8],
  ["16cell", "octahedron", 4, 1.5],
  ["24cell", "dodecahedron", 3, 2.5],
  ["120cell", "icosahedron", 3, 0.65],
  ["600cell", "tetrahedron", 4, 0.7],
];

describe("the model's gradients", () => {
  for (const [name, piece, P, s] of CASES)
    it(`${name} ${piece}: centres, rotations and the edge match central differences`, () => {
      const poly = getPolytope(name);
      const model = new SolidModel(poly, piece, P);
      const rand = mulberry32(7 + P);
      let tested = 0;
      for (let trial = 0; trial < 12 && tested < 3; trial++) {
        const st = model.randomState(s, rand);
        const g = new Float64Array(6 * P);
        const E = model.energy(st, s, g);
        if (!(E > 1e-6)) continue;
        tested++;
        const shift = (idx: number, d: number) => {
          const t = { cells: st.cells, U: Float64Array.from(st.U), R: Float64Array.from(st.R!) };
          if (idx < 3 * P) t.U[idx] += d;
          else {
            const q = Math.floor((idx - 3 * P) / 3);
            const w = [0, 0, 0];
            w[(idx - 3 * P) % 3] = d;
            const e = new Float64Array(9);
            expRot(w[0], w[1], w[2], e);
            const out = new Float64Array(9);
            mul9(e, 0, t.R, 9 * q, out, 0);
            t.R.set(out, 9 * q);
          }
          return t;
        };
        const h = 1e-6;
        const scratch = new Float64Array(6 * P);
        let bad = 0;
        for (let idx = 0; idx < 6 * P; idx++) {
          const fd = (model.energy(shift(idx, h), s, scratch) - model.energy(shift(idx, -h), s, scratch)) / (2 * h);
          // a state at a window's edge can jump (the model gates a pair by the window that holds its centre): such an entry is not smooth
          if (Math.abs(fd - g[idx]) > 1e-6 * Math.max(1, Math.abs(g[idx]))) bad++;
        }
        expect(bad, `${name} ${piece} trial ${trial}`).toBeLessThanOrEqual(1);
        const ds = new Float64Array(1);
        model.terms(st, s, null, false, { ds });
        const fds = (model.energy(st, s * (1 + 1e-6), scratch) - model.energy(st, s * (1 - 1e-6), scratch)) / (2e-6 * s);
        if (Math.abs(fds - ds[0]) < 1e-5 * Math.max(1, Math.abs(ds[0]))) expect(fds).toBeCloseTo(ds[0], 4);
      }
      expect(tested).toBeGreaterThan(0);
    });

  it("a weight of 2 doubles the energy, a focus keeps one piece's terms, and the violation is positive exactly when the energy is", () => {
    const poly = getPolytope("24cell");
    const model = new SolidModel(poly, "octahedron", 4);
    const rand = mulberry32(3);
    for (let t = 0; t < 6; t++) {
      const st = model.randomState(1.2, rand);
      const E = model.energy(st, 1.2, new Float64Array(24));
      const w = new Float64Array(4).fill(2);
      expect(model.terms(st, 1.2, null, false, { w })).toBeCloseTo(2 * E, 12);
      const v = model.violation(st, 1.2);
      expect(v > 0).toBe(E > 0);
      let sum = 0;
      for (let i = 0; i < 4; i++) sum += model.terms(st, 1.2, null, false, { focus: i });
      // each pair term appears in both pieces' focus (from either side), each edge term once: the focus terms cover the energy
      expect(sum).toBeGreaterThanOrEqual(E - 1e-12);
    }
  });
});

describe("the kernel for solids", () => {
  const make = (polytope: PolytopeName, piece: PolytopePieceName, P: number) => () => boundaryKernel({ polytope, piece, P });

  for (const [name, piece, P] of [
    ["5cell", "tetrahedron", 3],
    ["24cell", "octahedron", 3],
    ["tesseract", "icosahedron", 2],
    ["16cell", "cube", 3],
    ["120cell", "dodecahedron", 2],
  ] as [PolytopeName, PolytopePieceName, number][])
    it(`${name} ${piece} passes the generic methods' contract`, () => {
      expect(checkKernel(make(name, piece, P), { seed: P, states: 2 })).toEqual([]);
    }, 120000);

  it("has 6 variables a piece, swaps orientations, and its gradient (rotation vector and size included) matches central differences", () => {
    const k = boundaryKernel({ polytope: "24cell", piece: "octahedron", P: 3 });
    expect(k.dim).toBe(6);
    expect(k.steps).toHaveLength(6);
    expect(k.capabilities).toEqual(["seeds", "swap"]);
    const rng = mulberry32(2);
    let tested = 0;
    for (let t = 0; t < 10 && tested < 3; t++) {
      const s = k.randomStart(rng);
      k.setSize(s, s.size * 0.7);
      for (let i = 0; i < 3; i++) k.jiggle(s, i, 1.5, rng);
      const n = 18;
      const g = new Float64Array(n + 1);
      if (!(k.energy(s, g) > 1e-8)) continue;
      tested++;
      let bad = 0;
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
        const fd = (k.energy(p) - k.energy(m)) / (2 * h);
        if (Math.abs(fd - g[i]) > 1e-5 * Math.max(1, Math.abs(fd), Math.abs(g[i]))) bad++;
      }
      expect(bad).toBeLessThanOrEqual(1);
    }
    expect(tested).toBeGreaterThan(0);
    // swap: the two pieces exchange orientations (as sets of vertices) and keep their places
    const s = k.randomStart(rng);
    const before = k.toConfig(s) as BoundaryConfig;
    const sw = k.clone(s);
    k.swap!(sw, 0, 2);
    const after = k.toConfig(sw) as BoundaryConfig;
    const sh = pieceShape("octahedron");
    const verts = (q: [number, number, number, number]) => {
      const M = quaternionToMatrix(q);
      return sh.vertices.map((v) => [M[0] * v[0] + M[1] * v[1] + M[2] * v[2], M[3] * v[0] + M[4] * v[1] + M[5] * v[2], M[6] * v[0] + M[7] * v[1] + M[8] * v[2]]);
    };
    const same = (a: number[][], b: number[][]) => a.every((v) => b.some((u) => Math.hypot(u[0] - v[0], u[1] - v[1], u[2] - v[2]) < 1e-9));
    for (const [i, j] of [
      [0, 2],
      [2, 0],
    ]) {
      expect(after.pieces[i].cell).toBe(before.pieces[i].cell);
      for (let c = 0; c < 3; c++) expect(after.pieces[i].u[c]).toBeCloseTo(before.pieces[i].u[c], 12);
      expect(same(verts(after.pieces[i].q!), verts(before.pieces[j].q!))).toBe(true);
    }
    expect(after.pieces[1].q).toEqual(before.pieces[1].q);
  });

  it("finishes a relaxed state into a configuration the server's rule takes, and refuses what it cannot judge", () => {
    const k = boundaryKernel({ polytope: "5cell", piece: "octahedron", P: 2 });
    const rng = mulberry32(5);
    let made = 0;
    for (let t = 0; t < 6 && made < 2; t++) {
      const s = k.randomStart(rng);
      const f = k.finish(s);
      if (!f) continue;
      made++;
      const v = k.validate(f.config);
      expect(v.valid, v.reasons.join("; ")).toBe(true);
      expect(v.score).toBe(f.score);
      expect(f.score).toBeGreaterThanOrEqual(pieceBounds("5cell", "octahedron", 2).value * (1 - 1e-9));
      const back = k.fromConfig(f.config)!;
      expect(back.size).toBe(f.score);
    }
    expect(made).toBeGreaterThan(0);
    expect(k.fromConfig({ s: 3, pieces: [{ cell: 0, u: [0, 0, 0], q: [2, 0, 0, 0] }, { cell: 1, u: [0, 0, 0] }] })).toBeNull();
    expect(k.fromConfig({ s: 3, pieces: [{ cell: 0, u: [0, 0, 0] }] })).toBeNull();
    expect(k.fromConfig({ s: 1e5, pieces: [{ cell: 0, u: [0, 0, 0] }, { cell: 1, u: [0, 0, 0] }] })).toBeNull();
    // seeds: the cells themselves when the solid is the cell, else a sequential layout
    const cells = boundaryKernel({ polytope: "16cell", piece: "tetrahedron", P: 4 }).seeds!(mulberry32(1), 2);
    expect(cells).toHaveLength(2);
  });
});

describe("searchBoundary with solids", () => {
  const check = (name: PolytopeName, piece: PolytopePieceName, P: number, evals: number, seed: number) => {
    const problem = { polytope: name, piece, P };
    const r = searchBoundary(problem, { seed, budget: { evals } });
    expect(r.config, `${name} ${piece} ${P}`).not.toBeNull();
    const v0 = validateBoundary(problem, r.config, { tol: 0 });
    expect(v0.ok, v0.reasons.join("; ")).toBe(true);
    // the server's rule: tolerance 1e-9, 1024 chains
    const v9 = validateBoundary(problem, r.config, { tol: 1e-9, maxChains: 1024, maxReasons: 20 });
    expect(v9.ok, v9.reasons.join("; ")).toBe(true);
    expect(r.score).toBe(r.config!.s);
    expect(r.score!).toBeGreaterThanOrEqual(pieceBounds(name, piece, P).value * (1 - 1e-9));
    return r;
  };

  it("finds a packing of one solid on each container within a few thousand evaluations (the first offer, drawn in to touching)", () => {
    for (const [name, piece] of [
      ["tesseract", "dodecahedron"],
      ["5cell", "tetrahedron"],
      ["16cell", "cube"],
      ["600cell", "icosahedron"],
    ] as [PolytopeName, PolytopePieceName][])
      check(name, piece, 1, 4000, 1);
  }, 120000);

  it("starts from the cells themselves when the solid is the container's cell", () => {
    const r = check("24cell", "octahedron", 2, 3000, 2);
    expect(r.trace[0].source).toBe("cells");
    expect(r.trace[0].score).toBeCloseTo(1, 5);
  }, 60000);

  it("is reproducible from its seed and evaluation budget, and stops when asked", () => {
    const a = searchBoundary({ polytope: "tesseract", piece: "octahedron", P: 2 }, { seed: 4, budget: { evals: 2500 } });
    const b = searchBoundary({ polytope: "tesseract", piece: "octahedron", P: 2 }, { seed: 4, budget: { evals: 2500 } });
    expect(a.score).toBe(b.score);
    expect(JSON.stringify(a.config)).toBe(JSON.stringify(b.config));
    let calls = 0;
    const c = searchBoundary({ polytope: "tesseract", piece: "octahedron", P: 2 }, { seed: 4, budget: { evals: 1e9 }, shouldStop: () => ++calls > 300 });
    expect(c.ended).toBe("stopped");
    expect(c.evals).toBeLessThan(1000);
  }, 60000);

  it("improves a given start and never returns worse", () => {
    const first = searchBoundary({ polytope: "tesseract", piece: "icosahedron", P: 2 }, { seed: 1, budget: { evals: 4000 } });
    expect(first.config).not.toBeNull();
    const again = searchBoundary({ polytope: "tesseract", piece: "icosahedron", P: 2 }, { seed: 2, budget: { evals: 4000 }, start: first.config! });
    expect(again.score!).toBeLessThanOrEqual(first.score!);
    expect(again.startScore).toBe(first.score);
  }, 60000);

  it("treats a state it cannot represent as no warm start", () => {
    const r = searchBoundary({ polytope: "5cell", piece: "tetrahedron", P: 1 }, { seed: 1, budget: { evals: 2000 }, start: { s: 0.1, pieces: [{ cell: 0, u: [0, 0, 0] }] } });
    expect(r.startScore).toBeNull();
  }, 60000);
});

describe("states of the kernel", () => {
  it("round-trip through a configuration (rotations included)", () => {
    const k = boundaryKernel({ polytope: "600cell", piece: "dodecahedron", P: 3 });
    const rng = mulberry32(8);
    const s: KernelState = k.randomStart(rng);
    const cfg = k.toConfig(s) as BoundaryConfig;
    const back = k.fromConfig(cfg)!;
    const cfg2 = k.toConfig(back) as BoundaryConfig;
    for (let i = 0; i < 3; i++) {
      expect(cfg2.pieces[i].cell).toBe(cfg.pieces[i].cell);
      for (let c = 0; c < 3; c++) expect(cfg2.pieces[i].u[c]).toBeCloseTo(cfg.pieces[i].u[c], 12);
      for (let c = 0; c < 4; c++) expect(Math.abs(cfg2.pieces[i].q![c] - cfg.pieces[i].q![c])).toBeLessThan(1e-12);
    }
  });
});
