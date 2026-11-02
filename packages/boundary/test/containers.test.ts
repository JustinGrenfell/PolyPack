/**
 * The six containers of the boundary space (tesseract, 5-cell, 16-cell, 24-cell, 120-cell, 600-cell): the catalogue's polytopes
 * against the literature's counts and angles, their facts (bounds.ts), their ids, the lower bounds, and a gluing check on each.
 */
import { describe, expect, it } from "vitest";
import {
  CONTAINER_MAX_P,
  CONTAINER_PIECES,
  POLYTOPES,
  ballBounds,
  containerFacts,
  formatBoundaryId,
  getPolytope,
  isoApply,
  isoFromAffine,
  parseBoundaryId,
  toAmbient,
} from "../src/index";
import type { Polytope, PolytopeName } from "../src/index";

const DEG = 180 / Math.PI;
// vertices, edges, cells, ridges; cone angle of the report (section 2.5) and the cells around an edge
const LIT: Record<PolytopeName, { V: number; E: number; C: number; R: number; cone: number; around: number }> = {
  "5cell": { V: 5, E: 10, C: 5, R: 10, cone: 211.58633809652784, around: 3 },
  tesseract: { V: 16, E: 32, C: 8, R: 24, cone: 270, around: 3 },
  "16cell": { V: 8, E: 24, C: 16, R: 32, cone: 282.1151174620371, around: 4 },
  "24cell": { V: 24, E: 96, C: 24, R: 96, cone: 328.4136619034721, around: 3 },
  "120cell": { V: 600, E: 1200, C: 120, R: 720, cone: 349.6951535312339, around: 3 },
  "600cell": { V: 120, E: 720, C: 600, R: 1200, cone: 352.6438968275464, around: 5 },
};

describe("the containers", () => {
  it("lists the six, the tesseract first, and the pieces each takes", () => {
    expect([...POLYTOPES]).toEqual(["tesseract", "5cell", "16cell", "24cell", "120cell", "600cell"]);
    expect(CONTAINER_PIECES.tesseract).toEqual(["ball", "cube"]);
    for (const n of POLYTOPES.slice(1)) expect(CONTAINER_PIECES[n]).toEqual(["ball"]);
    expect(CONTAINER_MAX_P.tesseract).toBe(16);
    for (const n of POLYTOPES.slice(1)) expect(CONTAINER_MAX_P[n]).toBe(24);
  });

  for (const name of POLYTOPES) {
    const lit = LIT[name];
    describe(name, () => {
      const P: Polytope = getPolytope(name);

      it("has the literature's counts, unit edges and its centroid at the origin", () => {
        expect(P.vertices).toHaveLength(lit.V);
        expect(P.cells).toHaveLength(lit.C);
        expect(P.edges).toHaveLength(lit.E);
        expect(P.cells.reduce((n, c) => n + c.facets.length, 0)).toBe(2 * lit.R);
        for (const [a, b] of P.edges)
          expect(Math.hypot(...P.vertices[a].map((x, k) => x - P.vertices[b][k]))).toBeCloseTo(1, 12);
        for (let k = 0; k < 4; k++) expect(P.vertices.reduce((s, v) => s + v[k], 0) / lit.V).toBeCloseTo(0, 12);
        // no two vertices closer than an edge
        let least = Infinity;
        for (let i = 0; i < 5; i++)
          for (let j = i + 1; j < lit.V; j++)
            least = Math.min(least, Math.hypot(...P.vertices[i].map((x, k) => x - P.vertices[j][k])));
        expect(least).toBeGreaterThan(1 - 1e-12);
      });

      it("has the report's cone angle at every edge", () => {
        expect(P.hingeAngle * DEG).toBeCloseTo(lit.cone, 9);
        for (const a of P.edgeAngles) expect(a * DEG).toBeCloseTo(lit.cone, 9);
      });

      it("has gluings that are isometries fixing their ridge, in inverse pairs", () => {
        const s = 1.7;
        for (const c of P.cells.filter((_, i) => i % Math.max(1, Math.floor(lit.C / 12)) === 0))
          c.facets.forEach((f, fi) => {
            const A = f.toNeighbor.A;
            for (let i = 0; i < 3; i++)
              for (let j = 0; j < 3; j++)
                expect(A[0][i] * A[0][j] + A[1][i] * A[1][j] + A[2][i] * A[2][j]).toBeCloseTo(i === j ? 1 : 0, 12);
            const nb = P.cells[f.neighbor];
            expect(nb.facets[f.neighborFacet].neighbor).toBe(c.index);
            expect(nb.facets[f.neighborFacet].neighborFacet).toBe(fi);
            const g = isoFromAffine(f.toNeighbor, s);
            for (const vi of f.vertices) {
              const y2 = isoApply(
                g,
                c.vertices[vi].map((x) => x * s),
              );
              const x1 = toAmbient(
                P,
                c.index,
                c.vertices[vi].map((x) => x * s),
                s,
              );
              const x2 = toAmbient(P, nb.index, y2, s);
              for (let k = 0; k < 4; k++) expect(x2[k]).toBeCloseTo(x1[k], 12);
            }
          });
      });

      it("has facts that agree with its geometry", () => {
        const f = containerFacts(name);
        expect(f.cells).toBe(lit.C);
        expect(f.cellsAtEdge).toBe(lit.around);
        expect(f.coneAngle).toBeCloseTo(lit.cone, 9);
        // the cell volume by the pyramids over its facets
        let vol = 0;
        const c0 = P.cells[0];
        for (const fa of c0.facets) {
          const poly = fa.vertices.map((i) => c0.vertices[i]);
          let ax = 0,
            ay = 0,
            az = 0;
          for (let i = 0; i < poly.length; i++) {
            const a = poly[i],
              b = poly[(i + 1) % poly.length];
            ax += (a[1] - b[1]) * (a[2] + b[2]);
            ay += (a[2] - b[2]) * (a[0] + b[0]);
            az += (a[0] - b[0]) * (a[1] + b[1]);
          }
          vol += (fa.inside.h * Math.hypot(ax, ay, az)) / 6;
        }
        expect(vol).toBeCloseTo(f.cellVolume, 12);
        expect(f.boundaryVolume).toBeCloseTo(lit.C * vol, 10);
        // inradius: the facet planes' distance from the cell's centre; midradius: the centre's distance to the nearest edge
        for (const fa of c0.facets) expect(fa.inside.h).toBeCloseTo(f.inradius, 12);
        let near = Infinity;
        for (const [a, b] of c0.hinges) {
          const p = c0.vertices[a],
            q = c0.vertices[b];
          const e = q.map((x, k) => x - p[k]);
          const t = -(p[0] * e[0] + p[1] * e[1] + p[2] * e[2]) / (e[0] * e[0] + e[1] * e[1] + e[2] * e[2]);
          near = Math.min(near, Math.hypot(...p.map((x, k) => x + t * e[k])));
        }
        expect(near).toBeCloseTo(f.midradius, 12);
      });
    });
  }
});

describe("problem ids", () => {
  it("accepts balls on every container and cubes on the tesseract only", () => {
    for (const n of POLYTOPES) {
      const id = `boundary/${n}/ball/7`;
      expect(parseBoundaryId(id)).toEqual({ polytope: n, piece: "ball", P: 7 });
      expect(formatBoundaryId({ polytope: n, piece: "ball", P: 7 })).toBe(id);
    }
    expect(parseBoundaryId("boundary/tesseract/cube/12")).toEqual({ polytope: "tesseract", piece: "cube", P: 12 });
    for (const n of POLYTOPES.slice(1)) {
      expect(typeof parseBoundaryId(`boundary/${n}/cube/3`)).toBe("string");
      expect(() => formatBoundaryId({ polytope: n, piece: "cube", P: 3 })).toThrow();
    }
    expect(typeof parseBoundaryId("boundary/9cell/ball/3")).toBe("string");
    expect(typeof parseBoundaryId("boundary/24cell/ball/03")).toBe("string");
  });
});

describe("lower bounds of balls", () => {
  it("are the volume and the one-piece bounds, rounded down", () => {
    // one ball: 1/(2 rho) with rho the centre-to-edge distance of the cell
    const one: Record<string, number> = {
      "5cell": Math.SQRT2,
      "16cell": Math.SQRT2,
      "600cell": Math.SQRT2,
      "24cell": 1,
      "120cell": 2 / (3 + Math.sqrt(5)),
    };
    for (const [n, v] of Object.entries(one)) {
      const b = ballBounds(n as PolytopeName, 1);
      expect(b.onePiece).toBeLessThanOrEqual(v);
      expect(b.onePiece).toBeGreaterThan(v * (1 - 1e-14));
      expect(b.kind).toBe("one piece");
    }
    // the volume bound takes over with P: s^3 V = P pi / 6
    for (const n of POLYTOPES.slice(1)) {
      const f = containerFacts(n);
      const b = ballBounds(n, 24);
      expect(b.volume).toBeCloseTo(Math.cbrt((24 * Math.PI) / 6 / f.boundaryVolume), 12);
      expect(b.value).toBe(Math.max(b.volume, b.onePiece));
      expect(ballBounds(n, 24).volume).toBeGreaterThan(ballBounds(n, 12).volume);
    }
  });
});
