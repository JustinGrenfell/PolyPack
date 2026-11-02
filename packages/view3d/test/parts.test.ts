import { readFileSync } from "node:fs";
import { describe, expect, it } from "vitest";
import { cellPlanes, icosphere, partMesh, toPlane } from "../src/clip";
import { FIXTURES } from "../demo/fixtures";
import { validateBoundary } from "@polypack/boundary";
import { BALL_RADIUS, getPolytope, pieceParts, standardPiece } from "./engine";

/** Volume of a closed triangle soup (divergence theorem). */
function volume(pos: Float32Array): number {
  let v = 0;
  for (let i = 0; i < pos.length; i += 9) {
    const [ax, ay, az, bx, by, bz, cx, cy, cz] = pos.subarray(i, i + 9);
    v += (ax * (by * cz - bz * cy) - ay * (bx * cz - bz * cx) + az * (bx * cy - by * cx)) / 6;
  }
  return v;
}

const poly = getPolytope("tesseract");

describe("drawn parts of the engine's pieces", () => {
  // the parts of a piece partition it: their clipped meshes add up to the whole piece
  for (const fx of FIXTURES) {
    it(`add up to whole pieces: ${fx.label}`, () => {
      const shape = standardPiece(fx.piece, poly.d);
      const parts = pieceParts(poly, fx.piece, fx.config);
      const total = new Array<number>(fx.config.pieces.length).fill(0);
      for (const p of parts) {
        const planes = [...cellPlanes(poly.cells[p.cell], fx.config.s), ...p.clip.map((h) => toPlane(h))];
        const m = partMesh(shape, p.place, planes);
        if (m) total[p.piece] += volume(m.positions);
      }
      // A ball is drawn as the icosphere of detail 4, its cut faces as exact 96-gons: where they meet
      // the surface leaves gaps of 0.2% of the volume. One part missing or counted twice would be
      // 5% or more (a cap of the ball beyond a cell face at s = sqrt 2). Cubes are exact.
      // (the engine's ball parts carry the radius 1/2 in their place: the icosphere scaled by it)
      const whole = fx.piece === "ball" ? BALL_RADIUS ** 3 * volume(Float32Array.from(icosphere(4))) : 1;
      const tol = fx.piece === "ball" ? 5e-3 : 2e-4;
      for (const t of total) expect(Math.abs(t - whole) / whole).toBeLessThan(tol);
    });
  }
});

// The Platonic solids on the 4-polytopes (docs/spaces/BOUNDARY.md 4C): the drawing clips the standard solid by the cells' planes and the windows, as for
// the cube, so a solid's parts add up to the solid whatever its shape, rotation and container (the test packings of python/boundary/fixtures).
describe("drawn parts of the Platonic solids", () => {
  const packings = JSON.parse(readFileSync(new URL("../../../python/boundary/fixtures/platonic_packings.json", import.meta.url), "utf8")).packings as {
    id: string;
    polytope: string;
    piece: string;
    config: { s: number; pieces: { cell: number; u: number[]; q?: number[] }[] };
  }[];
  const VOLUME: Record<string, number> = {
    tetrahedron: Math.SQRT2 / 12,
    cube: 1,
    octahedron: Math.SQRT2 / 3,
    dodecahedron: (15 + 7 * Math.sqrt(5)) / 4,
    icosahedron: (5 * (3 + Math.sqrt(5))) / 12,
  };
  for (const r of packings.filter((p) => p.config.pieces.length <= 3))
    it(`add up to whole solids: ${r.id}`, () => {
      const pl = getPolytope(r.polytope as never);
      const shape = standardPiece(r.piece as never, pl.d);
      const parts = pieceParts(pl, r.piece as never, r.config as never);
      const total = new Array<number>(r.config.pieces.length).fill(0);
      for (const p of parts) {
        const planes = [...cellPlanes(pl.cells[p.cell], r.config.s), ...p.clip.map((h) => toPlane(h))];
        const m = partMesh(shape, p.place, planes);
        if (m) total[p.piece] += volume(m.positions);
      }
      for (const t of total) expect(Math.abs(t - VOLUME[r.piece]) / VOLUME[r.piece]).toBeLessThan(2e-4);
    });
});

// Random solids across ridges, in the engine's valid poses: the drawn parts add up to the solid (centres within a hair of a facet, so that a solid's face
// or vertex lies in or touches a cell's plane and the clipping meets the near-degenerate cases the gluings' rounding of 1e-10 makes).
describe("drawn parts of random Platonic solids across ridges", () => {
  const VOL: Record<string, number> = {
    tetrahedron: Math.SQRT2 / 12,
    cube: 1,
    octahedron: Math.SQRT2 / 3,
    dodecahedron: (15 + 7 * Math.sqrt(5)) / 4,
    icosahedron: (5 * (3 + Math.sqrt(5))) / 12,
  };
  let seed = 42;
  const rand = () => ((seed = (Math.imul(seed, 1664525) + 1013904223) >>> 0) / 4294967296);
  for (const name of ["tesseract", "5cell", "16cell", "24cell", "120cell", "600cell"] as const)
    for (const piece of ["tetrahedron", "cube", "octahedron", "dodecahedron", "icosahedron"] as const) {
      if (name === "tesseract" && piece === "cube") continue;
      it(`add up to whole solids: ${piece} on the ${name}`, () => {
        const pl = getPolytope(name);
        const shape = standardPiece(piece, pl.d);
        let tested = 0;
        for (let t = 0; t < 300 && tested < 3; t++) {
          const s = 3.5 * { tetrahedron: 0.612, cube: 0.866, octahedron: 0.707, dodecahedron: 1.4, icosahedron: 0.951 }[piece];
          const cell = Math.floor(rand() * pl.cells.length);
          const f = pl.cells[cell].facets[Math.floor(rand() * pl.cells[cell].facets.length)];
          const u = f.inside.n.map((x) => x * (f.inside.h * s - rand() * 0.3));
          const q = [rand() - 0.5, rand() - 0.5, rand() - 0.5, rand() - 0.5];
          const n = Math.hypot(...q);
          const config = { s, pieces: [{ cell, u, q: q.map((x) => x / n) as [number, number, number, number] }] };
          if (!validateBoundary({ polytope: name, piece, P: 1 }, config, { tol: 1e-9 }).ok) continue;
          const parts = pieceParts(pl, piece, config);
          if (parts.length < 2) continue;
          tested++;
          let total = 0;
          for (const p of parts) {
            const planes = [...cellPlanes(pl.cells[p.cell], s), ...p.clip.map((h) => toPlane(h))];
            const m = partMesh(shape, p.place, planes);
            if (m) total += volume(m.positions);
          }
          expect(Math.abs(total - VOL[piece]) / VOL[piece], `${name} ${piece} #${t}`).toBeLessThan(2e-4);
        }
        expect(tested).toBeGreaterThan(0);
      });
    }
});
