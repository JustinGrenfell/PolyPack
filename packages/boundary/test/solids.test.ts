/**
 * The Platonic solids as pieces (docs/spaces/BOUNDARY.md section 4C): the shapes, the engine's parts across a ridge, overlap and no-overlap at
 * known distances, the work caps (a check that cannot finish is refused), pose round trips, bounds, and the adversarial cases where the float
 * validator must never accept what the exact certifier (python/boundary/test_platonic.py) refuses.
 */
import { readFileSync } from "node:fs";
import { describe, expect, it } from "vitest";
import type { BoundaryConfig, PolytopeName, PolytopePieceName } from "../src/index";
import { MAX_PAIR_TESTS } from "../src/validate";
import {
  CONTAINER_MAX_P,
  POLYTOPE_PIECES,
  boundaryCheckMs,
  boundaryIsHeavy,
  boundaryVolume,
  containerFacts,
  getPolytope,
  isoFromAffine,
  matrixToQuaternion,
  parseBoundaryId,
  partMesh,
  pieceBounds,
  pieceParts,
  pieceShape,
  polyhedronVolume,
  quaternionToMatrix,
  rotationOnto,
  standardPiece,
  validateBoundary,
} from "../src/index";
import { polyhedronFrom, polyhedronVertices } from "../src/polyhedron";
import { SolidModel } from "../src/search/solids";
import { rng } from "./helpers";

const PHI = (1 + Math.sqrt(5)) / 2;
const SOLIDS = POLYTOPE_PIECES;
const NAMES: PolytopeName[] = ["tesseract", "5cell", "16cell", "24cell", "120cell", "600cell"];
/** The solid that is the cell of the container, if any. */
const CELL_OF: Partial<Record<PolytopeName, PolytopePieceName>> = { "5cell": "tetrahedron", "16cell": "tetrahedron", "600cell": "tetrahedron", "24cell": "octahedron", "120cell": "dodecahedron", tesseract: "cube" };

/** Each piece in the cell of the same shape, centred and turned onto it, at edge s. */
function cellsConfig(name: PolytopeName, piece: PolytopePieceName, s: number, count?: number): BoundaryConfig {
  const poly = getPolytope(name);
  const sh = pieceShape(piece);
  return {
    s,
    pieces: poly.cells.slice(0, count ?? poly.cells.length).map((c) => ({ cell: c.index, u: [0, 0, 0], q: matrixToQuaternion(rotationOnto(sh, c.vertices as [number, number, number][])!) })),
  };
}

describe("the five solids", () => {
  // vertices, edges, faces, volume, inradius, circumradius at unit edge (Euler: V - E + F = 2)
  const LIT: Record<PolytopePieceName, { V: number; E: number; F: number; vol: number; r: number; R: number; rho: number }> = {
    tetrahedron: { V: 4, E: 6, F: 4, vol: Math.SQRT2 / 12, r: Math.sqrt(6) / 12, R: Math.sqrt(6) / 4, rho: Math.SQRT2 / 4 },
    cube: { V: 8, E: 12, F: 6, vol: 1, r: 0.5, R: Math.sqrt(3) / 2, rho: Math.SQRT1_2 },
    octahedron: { V: 6, E: 12, F: 8, vol: Math.SQRT2 / 3, r: Math.sqrt(6) / 6, R: Math.SQRT1_2, rho: 0.5 },
    dodecahedron: { V: 20, E: 30, F: 12, vol: (15 + 7 * Math.sqrt(5)) / 4, r: Math.sqrt(250 + 110 * Math.sqrt(5)) / 20, R: (Math.sqrt(3) * PHI) / 2, rho: (PHI * PHI) / 2 },
    icosahedron: { V: 12, E: 30, F: 20, vol: (5 * (3 + Math.sqrt(5))) / 12, r: (PHI * PHI) / (2 * Math.sqrt(3)), R: Math.sqrt(PHI * Math.sqrt(5)) / 2, rho: PHI / 2 },
  };

  for (const name of SOLIDS)
    it(`${name} has the literature's counts, unit edges, a centred outward frame and its volume`, () => {
      const sh = pieceShape(name);
      const L = LIT[name];
      expect(sh.vertices).toHaveLength(L.V);
      expect(sh.edges).toHaveLength(L.E);
      expect(sh.faces).toHaveLength(L.F);
      for (const [a, b] of sh.edges) expect(Math.hypot(...sh.vertices[a].map((x, k) => x - sh.vertices[b][k]))).toBeCloseTo(1, 13);
      for (let k = 0; k < 3; k++) expect(sh.vertices.reduce((s, v) => s + v[k], 0)).toBeCloseTo(0, 13);
      for (const v of sh.vertices) expect(Math.hypot(...v)).toBeCloseTo(L.R, 13);
      sh.normals.forEach((n, f) => {
        for (const i of sh.faces[f]) expect(n[0] * sh.vertices[i][0] + n[1] * sh.vertices[i][1] + n[2] * sh.vertices[i][2]).toBeCloseTo(L.r, 13);
        for (const v of sh.vertices) expect(n[0] * v[0] + n[1] * v[1] + n[2] * v[2]).toBeLessThanOrEqual(L.r + 1e-13);
      });
      expect(sh.inradius).toBeCloseTo(L.r, 13);
      expect(sh.circumradius).toBeCloseTo(L.R, 13);
      expect(sh.midradius).toBeCloseTo(L.rho, 13);
      expect(sh.volume).toBeCloseTo(L.vol, 12);
      // the polyhedron built from the faces has that volume too (faces counter-clockwise from outside)
      expect(polyhedronVolume(polyhedronFrom(sh.vertices, sh.faces))).toBeCloseTo(L.vol, 12);
      // Euler
      expect(L.V - L.E + L.F).toBe(2);
      // the mirror maps the vertex set to itself and reverses handedness
      const M = sh.mirror;
      for (const v of sh.vertices) {
        const w = [M[0] * v[0] + M[1] * v[1] + M[2] * v[2], M[3] * v[0] + M[4] * v[1] + M[5] * v[2], M[6] * v[0] + M[7] * v[1] + M[8] * v[2]];
        expect(sh.vertices.some((u) => Math.hypot(u[0] - w[0], u[1] - w[1], u[2] - w[2]) < 1e-12)).toBe(true);
      }
      const det = M[0] * (M[4] * M[8] - M[5] * M[7]) - M[1] * (M[3] * M[8] - M[5] * M[6]) + M[2] * (M[3] * M[7] - M[4] * M[6]);
      expect(det).toBeCloseTo(-1, 12);
    });

  it("standardPiece gives each solid's vertices and faces (the cube as before)", () => {
    for (const name of SOLIDS) {
      const sp = standardPiece(name, 4);
      expect(sp.kind).toBe("polytope");
      if (sp.kind !== "polytope") continue;
      expect(sp.vertices).toHaveLength(LIT[name].V);
      expect(sp.faces).toHaveLength(LIT[name].F);
    }
    expect(() => standardPiece("sphere" as never, 4)).toThrow();
    expect(() => standardPiece("tetrahedron", 5)).toThrow();
  });

  it("rotationOnto turns each solid onto the cell of its container", () => {
    for (const name of NAMES) {
      const piece = CELL_OF[name]!;
      const sh = pieceShape(piece);
      for (const c of getPolytope(name).cells.slice(0, 6)) {
        const R = rotationOnto(sh, c.vertices as [number, number, number][])!;
        expect(R).not.toBeNull();
        for (const v of sh.vertices) {
          const w = [R[0] * v[0] + R[1] * v[1] + R[2] * v[2], R[3] * v[0] + R[4] * v[1] + R[5] * v[2], R[6] * v[0] + R[7] * v[1] + R[8] * v[2]];
          expect(c.vertices.some((u) => Math.hypot(u[0] - w[0], u[1] - w[1], u[2] - w[2]) < 1e-9)).toBe(true);
        }
      }
    }
  });
});

describe("the validator with solids", () => {
  it("accepts each container's cells, filled by the solid of their shape, at s = 1 (touching: valid at 1e-9, not strictly)", () => {
    for (const name of NAMES) {
      if (name === "tesseract") continue; // its cubes are the aligned ones
      const piece = CELL_OF[name]!;
      const count = name === "120cell" || name === "600cell" ? 24 : undefined;
      const P = count ?? getPolytope(name).cells.length;
      const ok = validateBoundary({ polytope: name, piece, P }, cellsConfig(name, piece, 1, count), { tol: 1e-9 });
      expect(ok.reasons, name).toEqual([]);
      expect(ok.ok).toBe(true);
      const more = validateBoundary({ polytope: name, piece, P }, cellsConfig(name, piece, 1 + 1e-5, count), { tol: 1e-9 });
      expect(more.ok, name).toBe(true);
      const less = validateBoundary({ polytope: name, piece, P }, cellsConfig(name, piece, 1 - 1e-5, count), { tol: 1e-9 });
      expect(less.ok, name).toBe(false);
      expect(less.reasons.join(" ")).toMatch(/total volume|passing through its interior|overlap/);
    }
  });

  it("sees an overlap in the neighbouring cell at known depths and none at the contact", () => {
    // two octahedra in adjacent cells of the 24-cell, each filling its cell; moving one towards the other by d overlaps them by d
    const poly = getPolytope("24cell");
    const base = cellsConfig("24cell", "octahedron", 1, 24);
    const a = poly.cells[0];
    const f = a.facets[0];
    const nb = f.neighbor;
    // the facet's inward normal in cell a's frame; its inside.n points out, so moving along it goes towards the neighbour
    for (const [d, ok] of [
      [0, true],
      [1e-12, true],
      [1e-6, false],
      [1e-3, false],
    ] as const) {
      const cfg: BoundaryConfig = { s: 1, pieces: base.pieces.map((p) => ({ ...p, u: p.u.slice() })) };
      const only = [cfg.pieces[0], cfg.pieces[nb]];
      only[0].u = f.inside.n.map((x) => x * d);
      const sub: BoundaryConfig = { s: 1, pieces: only };
      const v = validateBoundary({ polytope: "24cell", piece: "octahedron", P: 2 }, sub, { tol: 1e-9 });
      expect(v.ok, `d = ${d}: ${v.reasons[0]}`).toBe(ok);
      if (!ok) expect(v.reasons.join(" ")).toMatch(/overlap|passing through/);
    }
  });

  it("splits a solid across a ridge into parts that add up to its volume, in the right cells", () => {
    const rand = rng(11);
    for (const name of ["16cell", "24cell", "120cell"] as PolytopeName[])
      for (const piece of SOLIDS) {
        const poly = getPolytope(name);
        const sh = pieceShape(piece);
        let tested = 0;
        for (let t = 0; t < 200 && tested < 2; t++) {
          const cell = Math.floor(rand() * poly.cells.length);
          // a centre near a facet of its cell, a random turn, an edge large enough that a solid fits
          const s = 3.2 * sh.circumradius / (name === "120cell" ? 1 : 1);
          const f = poly.cells[cell].facets[Math.floor(rand() * poly.cells[cell].facets.length)];
          const u = f.inside.n.map((x) => x * (f.inside.h * s - 0.45 * rand() * sh.inradius));
          const q = matrixToQuaternion(quaternionToMatrix([rand() - 0.5, rand() - 0.5, rand() - 0.5, rand() - 0.5]));
          const config: BoundaryConfig = { s, pieces: [{ cell, u, q }] };
          const v = validateBoundary({ polytope: name, piece, P: 1 }, config, { tol: 1e-9 });
          if (!v.ok) continue;
          const parts = pieceParts(poly, piece, config);
          let total = 0;
          for (const part of parts) {
            const mesh = partMesh(poly, piece, part, s)!;
            total += polyhedronVolume(polyhedronFrom(mesh.vertices as [number, number, number][], mesh.faces));
            // inside its cell: every vertex satisfies the cell's facets
            for (const vtx of mesh.vertices) for (const fa of poly.cells[part.cell].facets) expect(fa.inside.n[0] * vtx[0] + fa.inside.n[1] * vtx[1] + fa.inside.n[2] * vtx[2]).toBeLessThanOrEqual(fa.inside.h * s + 1e-9);
          }
          expect(total).toBeCloseTo(sh.volume, 9);
          if (parts.length > 1) tested++;
          // the parts' placements are isometries
          for (const part of parts) {
            const T = isoFromAffine(part.place);
            const det = T.m[0] * (T.m[4] * T.m[8] - T.m[5] * T.m[7]) - T.m[1] * (T.m[3] * T.m[8] - T.m[5] * T.m[6]) + T.m[2] * (T.m[3] * T.m[7] - T.m[4] * T.m[6]);
            expect(Math.abs(Math.abs(det) - 1)).toBeLessThan(1e-12);
          }
        }
        expect(tested, `${name} ${piece}`).toBeGreaterThan(0);
      }
  });

  it("refuses an edge through a solid, at every depth beyond the tolerance", () => {
    // the centre of the 5-cell's tetrahedron is rho = sqrt2/4 = 0.3536 from its edges (unit edge): a solid whose inradius is more than that cannot sit there
    const poly = getPolytope("5cell");
    for (const piece of SOLIDS) {
      const sh = pieceShape(piece);
      const rho = 0.35355339059327373;
      const s = 1.2;
      const q = matrixToQuaternion(quaternionToMatrix([0.3, 0.5, -0.2, 0.7]));
      const centred: BoundaryConfig = { s, pieces: [{ cell: 0, u: [0, 0, 0], q }] };
      const v = validateBoundary({ polytope: "5cell", piece, P: 1 }, centred, { tol: 1e-9 });
      if (sh.inradius > rho * s + 1e-9) {
        expect(v.ok, piece).toBe(false);
        expect(v.reasons.join(" ")).toMatch(/edge of the 5cell passing through|total volume/);
      }
      expect(poly.cells[0].hinges.length).toBe(6);
    }
  });

  it("is refused, never accepted, when a work cap is hit", () => {
    const cfg = cellsConfig("16cell", "tetrahedron", 1 + 1e-5);
    const P = 16;
    expect(validateBoundary({ polytope: "16cell", piece: "tetrahedron", P }, cfg, { tol: 1e-9 }).ok).toBe(true);
    const chains = validateBoundary({ polytope: "16cell", piece: "tetrahedron", P }, cfg, { tol: 1e-9, maxChains: 3 });
    expect(chains.ok).toBe(false);
    expect(chains.reasons[0]).toMatch(/more than 3 chains of cells/);
    // two octahedra in one cell of the 24-cell at s = 2.5, face to face 0.9 apart (valid): the parts share a cell, so a pair is tested
    const n3 = 0.45 / Math.sqrt(3);
    const pairCfg: BoundaryConfig = { s: 2.5, pieces: [{ cell: 0, u: [-n3, -n3, -n3] }, { cell: 0, u: [n3, n3, n3] }] };
    const free = validateBoundary({ polytope: "24cell", piece: "octahedron", P: 2 }, pairCfg, { tol: 1e-9 });
    expect(free.ok, free.reasons[0]).toBe(true);
    const pairs = validateBoundary({ polytope: "24cell", piece: "octahedron", P: 2 }, pairCfg, { tol: 1e-9, maxPairs: 0 });
    expect(pairs.ok).toBe(false);
    expect(pairs.reasons[0]).toMatch(/more than 0 tests of pairs of parts/);
    // the default cap (MAX_PAIR_TESTS tests, about 0.2 s at the slowest measured cost) refuses a crowd of copies in one cell, with its reason
    expect(MAX_PAIR_TESTS).toBeLessThanOrEqual(20000);
    const copies = 250; // 250 x 249 / 2 = 31,125 overlapping pairs of bounding boxes
    const heap: BoundaryConfig = { s: 5, pieces: Array.from({ length: copies }, () => ({ cell: 0, u: [0, 0, 0] })) };
    const crowd = validateBoundary({ polytope: "24cell", piece: "octahedron", P: copies }, heap, { tol: 1e-9, maxReasons: 1e9 });
    expect(crowd.ok).toBe(false);
    expect(crowd.reasons.join(" ")).toMatch(new RegExp(`more than ${MAX_PAIR_TESTS} tests of pairs of parts`));
    // the cap on reasons never turns a failure into a pass
    const crowded: BoundaryConfig = { s: 1 + 1e-5, pieces: Array.from({ length: 16 }, () => cfg.pieces[0]) };
    const r = validateBoundary({ polytope: "16cell", piece: "tetrahedron", P }, crowded, { tol: 1e-9, maxReasons: 1 });
    expect(r.ok).toBe(false);
    expect(r.reasons).toHaveLength(1);
  });

  it("the tesseract's cube is the axis-aligned problem: q is absent or one of the cube's symmetries, a turned cube is refused", () => {
    const prob = { polytope: "tesseract", piece: "cube", P: 1 } as const;
    const at = (q?: number[]): BoundaryConfig => ({ s: 3, pieces: [q ? { cell: 0, u: [0, 0, 0], q: q as [number, number, number, number] } : { cell: 0, u: [0, 0, 0] }] });
    const h = Math.SQRT1_2;
    for (const q of [undefined, [1, 0, 0, 0], [-1, 0, 0, 0], [h, h, 0, 0], [0.5, 0.5, 0.5, 0.5], [0, h, h, 0]]) {
      const v = validateBoundary(prob, at(q), { tol: 1e-9 });
      expect(v.ok, JSON.stringify(q)).toBe(true);
    }
    // 45 degrees about z, and a turn about no axis: cubes in space, but not the aligned problem
    const a = Math.PI / 8;
    const turned = [Math.cos(a), 0, 0, Math.sin(a)] as [number, number, number, number];
    const skew = [0.6, 0.48, 0.64, 0] as [number, number, number, number];
    for (const q of [turned, skew, [h, h + 1e-6, 0, 0].map((x) => x / Math.hypot(h, h + 1e-6))]) {
      const v = validateBoundary(prob, at(q), { tol: 1e-9 });
      expect(v.ok, JSON.stringify(q)).toBe(false);
      expect(v.reasons[0]).toMatch(/axis-aligned/);
    }
    // the other solids on the tesseract, and the cube on the other containers, still take any rotation
    expect(validateBoundary({ polytope: "tesseract", piece: "octahedron", P: 1 }, at(turned), { tol: 1e-9 }).ok).toBe(true);
    expect(validateBoundary({ polytope: "24cell", piece: "cube", P: 1 }, at(turned), { tol: 1e-9 }).ok).toBe(true);
  });

  it("rejects a malformed solid: a bad quaternion, a rotation on a ball, wrong P, NaN", () => {
    const cfg = cellsConfig("5cell", "tetrahedron", 1.1, 2);
    const ok = validateBoundary({ polytope: "5cell", piece: "tetrahedron", P: 2 }, cfg).ok;
    expect(ok).toBe(true);
    const bad = (c: unknown, piece: PolytopePieceName | "ball" = "tetrahedron", P = 2) => validateBoundary({ polytope: "5cell", piece, P }, c);
    expect(bad({ s: 1.1, pieces: [{ cell: 0, u: [0, 0, 0], q: [2, 0, 0, 0] }, cfg.pieces[1]] }).ok).toBe(false);
    expect(bad({ s: 1.1, pieces: [{ cell: 0, u: [0, 0, 0], q: [NaN, 0, 0, 0] }, cfg.pieces[1]] }).ok).toBe(false);
    expect(bad(cfg, "ball").ok).toBe(false);
    expect(bad(cfg, "tetrahedron", 3).ok).toBe(false);
    expect(bad({ s: 1.1, pieces: [{ cell: 9, u: [0, 0, 0] }, cfg.pieces[1]] }).ok).toBe(false);
    // a solid on a d-cube is no problem
    expect(validateBoundary({ polytope: "6cube", piece: "tetrahedron", P: 1 }, { s: 3, pieces: [{ cell: 0, u: [0, 0, 0, 0, 0] }] }).reasons[0]).toMatch(/balls and cubes only/);
  });
});

describe("bounds", () => {
  it("are the volume and one-piece bounds, each below a packing that exists", () => {
    for (const name of NAMES)
      for (const piece of SOLIDS) {
        if (name === "tesseract" && piece === "cube") continue;
        const sh = pieceShape(piece);
        const facts = { boundaryVolume: boundaryVolume(getPolytope(name)) };
        for (const P of [1, 5, 24]) {
          const b = pieceBounds(name, piece, P);
          expect(b.value).toBe(Math.max(b.volume, b.onePiece));
          expect(b.volume).toBeCloseTo(Math.cbrt((P * sh.volume) / facts.boundaryVolume), 12);
          expect(b.volume).toBeLessThan(Math.cbrt((P * sh.volume) / facts.boundaryVolume) * (1 + 1e-14));
          // the one-piece bound is r / rho, rho the midradius of the container's cell, rounded down
          expect(b.onePiece).toBeCloseTo(sh.inradius / containerFacts(name).midradius, 13);
          expect(b.onePiece).toBeLessThanOrEqual(sh.inradius / containerFacts(name).midradius);
        }
      }
    // the cells themselves are packings at s = 1: every bound is below 1 there
    for (const name of ["5cell", "16cell", "24cell", "120cell", "600cell"] as PolytopeName[]) {
      const piece = CELL_OF[name]!;
      const poly = getPolytope(name);
      expect(pieceBounds(name, piece, poly.cells.length > 24 ? 24 : poly.cells.length).value).toBeLessThan(1);
    }
  });

  it("agree with the closed forms: the 5-cell's tetrahedron, the 120-cell's dodecahedron", () => {
    // 5-cell: rho = sqrt2/4, r = sqrt6/12  ->  r / rho = 1 / sqrt 3
    expect(pieceBounds("5cell", "tetrahedron", 1).onePiece).toBeCloseTo(1 / Math.sqrt(3), 14);
    // the 5 tetrahedra fill the 5-cell exactly at s = 1: the volume bound at P = 5 is 1
    expect(pieceBounds("5cell", "tetrahedron", 5).volume).toBeCloseTo(1, 14);
    expect(pieceBounds("16cell", "tetrahedron", 16).volume).toBeCloseTo(1, 14);
    expect(pieceBounds("24cell", "octahedron", 24).volume).toBeCloseTo(1, 14);
    // 120-cell: r = 1.1135, rho = phi^2 / 2
    expect(pieceBounds("120cell", "dodecahedron", 1).onePiece).toBeCloseTo(pieceShape("dodecahedron").inradius / ((PHI * PHI) / 2), 13);
    // the tesseract: r / (sqrt 2 / 2); a cube would be 1 / sqrt 2 as in section 5
    expect(pieceBounds("tesseract", "cube", 1).onePiece).toBeCloseTo(Math.SQRT1_2, 14);
  });

  it("price a solid at the measured per-piece times, heavy from the stated P", () => {
    expect(boundaryCheckMs({ polytope: "120cell", piece: "dodecahedron", P: 5 })).toBe(55);
    expect(boundaryIsHeavy({ polytope: "120cell", piece: "dodecahedron", P: 5 })).toBe(true);
    expect(boundaryIsHeavy({ polytope: "120cell", piece: "dodecahedron", P: 4 })).toBe(false);
    expect(boundaryIsHeavy({ polytope: "600cell", piece: "icosahedron", P: 24 })).toBe(false);
    expect(boundaryIsHeavy({ polytope: "5cell", piece: "cube", P: 9 })).toBe(true);
    expect(boundaryIsHeavy({ polytope: "5cell", piece: "cube", P: 8 })).toBe(false);
    // the tesseract's aligned cubes keep their price, balls theirs
    expect(boundaryCheckMs({ polytope: "tesseract", piece: "cube", P: 16 })).toBeCloseTo(4.8, 12);
    expect(boundaryCheckMs({ polytope: "24cell", piece: "ball", P: 24 })).toBeCloseTo(55.2, 12);
    expect(CONTAINER_MAX_P.tesseract).toBe(16);
  });
});

describe("orientations of the search's model", () => {
  it("round-trip through the configuration, and an improper orientation is the same piece turned by the mirror", () => {
    const rand = rng(5);
    for (const piece of SOLIDS) {
      const poly = getPolytope("16cell");
      const model = new SolidModel(poly, piece, 2);
      const sh = pieceShape(piece);
      const st = model.randomState(model.startEdge(), rand);
      const cfg = { s: model.startEdge(), pieces: model.pieces(st, model.startEdge()) };
      const back = model.stateOf(cfg)!;
      for (let k = 0; k < 18; k++) expect(back.R![k]).toBeCloseTo(st.R![k], 12);
      // an improper orientation
      const R = Array.from(st.R!.subarray(0, 9));
      const flipped = R.map((x, i) => (i >= 6 ? -x : x));
      const st2 = { cells: st.cells, U: st.U, R: Float64Array.from([...flipped, ...R]) };
      const out = model.pieces(st2, model.startEdge());
      const Rq = quaternionToMatrix(out[0].q!);
      const verts = (M: number[]) => sh.vertices.map((v) => [M[0] * v[0] + M[1] * v[1] + M[2] * v[2], M[3] * v[0] + M[4] * v[1] + M[5] * v[2], M[6] * v[0] + M[7] * v[1] + M[8] * v[2]]);
      const want = verts(flipped);
      for (const w of verts(Rq)) expect(want.some((u) => Math.hypot(u[0] - w[0], u[1] - w[1], u[2] - w[2]) < 1e-9)).toBe(true);
    }
  });
});

describe("the cross-check fixture", () => {
  it("still holds: the engine reproduces every verdict recorded for the exact certifier's test, and every lower bound", () => {
    const data = JSON.parse(readFileSync(new URL("../../../python/boundary/fixtures/platonic_crosscheck.json", import.meta.url), "utf8")) as {
      entries: { polytope: PolytopeName; piece: PolytopePieceName; P: number; tol: number; cases: { sign: number; config: BoundaryConfig; engine: boolean }[] }[];
      bounds: Record<string, number>;
    };
    expect(data.entries.length).toBeGreaterThan(50);
    for (const e of data.entries)
      for (const c of e.cases) {
        const v = validateBoundary({ polytope: e.polytope, piece: e.piece, P: e.P }, c.config, { tol: e.tol, maxChains: 4096 });
        expect(v.ok, `${e.polytope} ${e.piece} ${c.sign}: ${v.reasons[0]}`).toBe(c.engine);
        expect(c.engine).toBe(c.sign > 0);
      }
    expect(Object.keys(data.bounds)).toHaveLength(664);
    for (const [id, value] of Object.entries(data.bounds)) {
      const p = parseBoundaryId(id) as { polytope: PolytopeName; piece: PolytopePieceName; P: number };
      expect(pieceBounds(p.polytope, p.piece, p.P).value, id).toBe(value);
    }
  });
});
