/**
 * Balls on the boundary of the d-simplex, d-orthoplex and d-demicube, d = 5..10 (docs/spaces/BOUNDARY.md section 4C): the containers'
 * combinatorics, the development maps, the half-chain engine against random walks, known distances, the bounds, the float validator, the server
 * rule, the kernel contract and the seeds. The exact certifier (python/boundary/ndpoly.py) has its own tests and reads the engine's
 * numbers written by scripts/gen-nd-crosscheck.ts; here the engine reads the certifier's test packings (fixtures/nd_packings.json).
 */
import { readFileSync } from "node:fs";
import { fileURLToPath } from "node:url";
import { checkKernel } from "@polypack/methods/contract";
import { mulberry32 } from "@polypack/methods";
import { describe, expect, it } from "vitest";
import {
  BOUNDARY_HEAVY_MS,
  CONTAINER_MAX_P,
  CONTAINERS,
  ND_MAX_P,
  ND_NAMES,
  NDLimitError,
  boundaryDomain,
  boundaryIsHeavy,
  boundaryKernel,
  getND,
  isContainerName,
  measureND,
  ndBallVolume,
  ndBound,
  ndBoundKind,
  ndBoundaryVolume,
  ndCellCount,
  ndCentreEdge,
  ndCounts,
  ndDevelop,
  ndDirectionAtAnchor,
  ndExactLayout,
  ndHingeDistance,
  ndJoinDistance,
  ndPairEdge,
  ndParse,
  ndSequentialLayout,
  ndServerRule,
  newNDBudget,
  parseBoundaryId,
  searchBoundary,
  validateBoundary,
  validateND,
  walkND,
} from "../src/index";
import type { NDName } from "../src/index";

type Cfg = { s: number; pieces: { cell: number; u: number[] }[] };
const fixtures = (name: string) =>
  JSON.parse(
    readFileSync(fileURLToPath(new URL(`../../../python/boundary/fixtures/${name}`, import.meta.url)), "utf8"),
  );
const packings: { id: string; polytope: NDName; P: number; source: string; config: Cfg }[] =
  fixtures("nd_packings.json").packings;
const problem = (polytope: NDName, P: number) => ({ polytope, piece: "ball" as const, P });
const dot = (a: ArrayLike<number>, b: ArrayLike<number>) => {
  let s = 0;
  for (let i = 0; i < a.length; i++) s += a[i] * b[i];
  return s;
};
const zeros = (n: number) => new Array<number>(n).fill(0);
const mulberry = (seed: number) => mulberry32(seed);
const gauss = (rand: () => number) => Math.sqrt(-2 * Math.log(1 - rand())) * Math.cos(2 * Math.PI * rand());

/** A random point of a cell, inside it (shrunk toward the centre until it is). */
function randomPoint(name: NDName, c: number, s: number, rand: () => number, spread: number): number[] {
  const C = getND(name);
  const cl = C.cell(c);
  const u0 = Array.from({ length: C.n }, () => (gauss(rand) * spread * s) / Math.sqrt(C.n));
  for (let lam = 1, t = 0; t < 40; t++, lam *= 0.8) {
    const u = u0.map((v) => v * lam);
    const x = C.toAmbient(c, u, s);
    let ok = true;
    for (let k = 0; k < cl.nbr.length && ok; k++) {
      let v = -s * cl.hn[k];
      for (let r = 0; r < C.D; r++) v += cl.nn[k * C.D + r] * x[r];
      if (v > 0) ok = false;
    }
    if (ok) return u;
  }
  return zeros(C.n);
}
function randomTangent(name: NDName, c: number, len: number, rand: () => number): Float64Array {
  const C = getND(name);
  const dv = new Float64Array(C.D);
  for (let k = 0; k < C.n; k++) {
    const g = gauss(rand);
    for (let j = 0; j < C.D; j++) dv[j] += g * C.cell(c).frame[k][j];
  }
  const nn = Math.hypot(...dv);
  return dv.map((v) => (v * len) / nn);
}

describe("the containers", () => {
  it("have the facet counts of the literature [EST: Coxeter, Regular Polytopes]", () => {
    for (const name of ND_NAMES) {
      const p = ndParse(name)!;
      const expected =
        p.family === "simplex" ? p.d + 1 : p.family === "orthoplex" ? 2 ** p.d : 2 * p.d + 2 ** (p.d - 1);
      expect(getND(name).C, name).toBe(expected);
      expect(ndCellCount(p.family, p.d)).toBe(expected);
    }
    // the 5-demicube: 16 5-cells and 10 16-cells; its f-vector ends 120 3-faces, 26 4-faces
    expect(ndCounts("demicube", 5)).toEqual({ facets: 26, ridges: 120, simplexFacets: 16, otherFacets: 10 });
    expect(ndCounts("orthoplex", 10).facets).toBe(1024);
    expect(ndCounts("demicube", 10).facets).toBe(532);
    expect(ndCounts("demicube", 10)).toMatchObject({ simplexFacets: 512, otherFacets: 20 });
    expect(ndCounts("simplex", 7)).toMatchObject({ facets: 8, ridges: 28 });
    expect(ndCounts("orthoplex", 5)).toMatchObject({ facets: 32, ridges: 80 });
  });

  it("are convex polytopes of unit edge whose cells, ridges and hinges match a brute-force reading of the vertex set (d = 5, 6, 7)", () => {
    for (const name of ND_NAMES.filter((n) => /^[567]/.test(n))) {
      const C = getND(name);
      const V = C.vertices(1);
      // unit edge: the closest pair of vertices is 1 apart
      let minD = Infinity;
      for (let i = 0; i < V.length; i++)
        for (let j = i + 1; j < V.length; j++) minD = Math.min(minD, Math.hypot(...V[i].map((x, r) => x - V[j][r])));
      expect(minD, name).toBeCloseTo(1, 12);
      const cells =
        C.C <= 70
          ? Array.from({ length: C.C }, (_, i) => i)
          : [0, 1, 2, 5, 7, 11, 13, 20, 33, 40, 41, 63].filter((c) => c < C.C);
      const tight = (c: number) => {
        const cl = C.cell(c);
        return V.map((v, i) => [Math.abs(dot(cl.n, v) - cl.h) < 1e-9, i] as const)
          .filter((t) => t[0])
          .map((t) => t[1]);
      };
      const rank = (ids: number[]): number => {
        if (ids.length < 2) return 0;
        const rows = ids.slice(1).map((i) => Array.from(V[i], (x, r) => x - V[ids[0]][r]));
        let rk = 0;
        for (let col = 0; col < C.D && rk < rows.length; col++) {
          let piv = rk;
          for (let i = rk + 1; i < rows.length; i++) if (Math.abs(rows[i][col]) > Math.abs(rows[piv][col])) piv = i;
          if (Math.abs(rows[piv][col]) < 1e-9) continue;
          [rows[rk], rows[piv]] = [rows[piv], rows[rk]];
          for (let i = rk + 1; i < rows.length; i++) {
            const f = rows[i][col] / rows[rk][col];
            for (let r = col; r < C.D; r++) rows[i][r] -= f * rows[rk][r];
          }
          rk++;
        }
        return rk;
      };
      for (const c of cells) {
        const cl = C.cell(c);
        // every vertex is inside; the cell's vertices span n dimensions
        for (const v of V) expect(dot(cl.n, v) - cl.h).toBeLessThan(1e-9);
        const tv = tight(c);
        expect(rank(tv), `${name} cell ${c}`).toBe(C.n);
        // the neighbours are exactly the cells sharing an (n-1)-dimensional face
        const nbrBrute: number[] = [];
        for (let g = 0; g < C.C; g++) {
          if (g === c) continue;
          const common = tv.filter((i) => Math.abs(dot(C.cell(g).n, V[i]) - C.cell(g).h) < 1e-9);
          if (common.length >= C.n && rank(common) === C.n - 1) nbrBrute.push(g);
        }
        expect(
          [...cl.nbr].sort((a, b) => a - b),
          `${name} cell ${c} neighbours`,
        ).toEqual(nbrBrute);
        // the ridge points lie on both planes
        for (let k = 0; k < cl.nbr.length; k++) {
          const q = cl.rp.subarray(k * C.D, (k + 1) * C.D);
          expect(Math.abs(dot(cl.n, q) - cl.h)).toBeLessThan(1e-12);
          expect(Math.abs(dot(C.cell(cl.nbr[k]).n, q) - C.cell(cl.nbr[k]).h)).toBeLessThan(1e-12);
        }
        // the hinge pairs are exactly the pairs of ridges meeting in an (n-2)-dimensional face
        if (C.C <= 70 || c < 3) {
          const pairs = new Set<string>();
          for (let p = 0; p < cl.hinges.length; p += 2) pairs.add(`${cl.hinges[p]},${cl.hinges[p + 1]}`);
          const brute = new Set<string>();
          for (let k = 0; k < cl.nbr.length; k++)
            for (let l = k + 1; l < cl.nbr.length; l++) {
              const common = tv.filter(
                (i) =>
                  Math.abs(dot(C.cell(cl.nbr[k]).n, V[i]) - C.cell(cl.nbr[k]).h) < 1e-9 &&
                  Math.abs(dot(C.cell(cl.nbr[l]).n, V[i]) - C.cell(cl.nbr[l]).h) < 1e-9,
              );
              if (common.length >= C.n - 1 && rank(common) === C.n - 2) brute.add(`${k},${l}`);
            }
          expect(pairs, `${name} cell ${c} hinges`).toEqual(brute);
        }
      }
    }
  }, 300_000);

  it("meet at the dihedral angles of the literature: cos = 1/d (simplex), (d-2)/d (orthoplex); the demicube's facets are square to each other and cross a simplex facet at cos = 1/sqrt d", () => {
    for (const name of ND_NAMES) {
      const p = ndParse(name)!;
      const C = getND(name);
      const cosines = new Set<string>();
      for (let c = 0; c < Math.min(C.C, 20); c++) {
        const cl = C.cell(c);
        // the interior dihedral angle between cells c and g has cosine -n_c . n_g
        for (let k = 0; k < cl.nbr.length; k++)
          cosines.add(Math.abs(dot(cl.n, cl.nn.subarray(k * C.D, (k + 1) * C.D))).toFixed(10));
      }
      if (p.family === "simplex") expect([...cosines], name).toEqual([(1 / p.d).toFixed(10)]);
      if (p.family === "orthoplex") expect([...cosines], name).toEqual([((p.d - 2) / p.d).toFixed(10)]);
      if (p.family === "demicube")
        expect([...cosines].sort(), name).toEqual([(0).toFixed(10), (1 / Math.sqrt(p.d)).toFixed(10)].sort());
    }
  });
});

describe("the development maps", () => {
  it("are isometries that fix the ridge and lay the neighbour on the far side (the reflection with normal n1 + n2), and agree with the rotation of the two normals", () => {
    const rand = mulberry(11);
    for (const name of ND_NAMES) {
      const C = getND(name);
      const D = C.D;
      for (let t = 0; t < 6; t++) {
        const c = Math.floor(rand() * C.C);
        const cl = C.cell(c);
        const k = Math.floor(rand() * cl.nbr.length);
        const g = C.cell(cl.nbr[k]);
        const m = Array.from(cl.n, (v, j) => v + g.n[j]);
        const mm = dot(m, m);
        const q = cl.rp.subarray(k * D, (k + 1) * D);
        const sigma = (x: ArrayLike<number>) => {
          const f =
            (2 *
              dot(
                m,
                Array.from(x, (v, j) => v - q[j]),
              )) /
            mm;
          return Array.from(x, (v, j) => v - f * m[j]);
        };
        // involution, isometry
        const a = Array.from({ length: D }, () => gauss(rand));
        const b = Array.from({ length: D }, () => gauss(rand));
        expect(Math.hypot(...sigma(a).map((v, j) => v - sigma(b)[j]))).toBeCloseTo(
          Math.hypot(...a.map((v, j) => v - b[j])),
          10,
        );
        sigma(sigma(a)).forEach((v, j) => expect(v).toBeCloseTo(a[j], 10));
        // fixes the ridge point; maps the plane of g onto the plane of c
        sigma(q).forEach((v, j) => expect(v).toBeCloseTo(q[j], 12));
        const og = Array.from(g.o);
        const img = sigma(og);
        expect(Math.abs(dot(cl.n, img) - cl.h)).toBeLessThan(1e-12);
        // the far side: the image of g's centre is beyond the ridge, n_g . x > h_g, by exactly h_c - n_c . o_g ... = the height of g's centre above its own ridge plane
        expect(dot(g.n, img) - g.h).toBeGreaterThan(1e-9);
        // the other reflection through the ridge (normal n1 - n2) folds the neighbour back over c: its image is on c's own side
        const m2 = Array.from(cl.n, (v, j) => v - g.n[j]);
        const f2 =
          (2 *
            dot(
              m2,
              og.map((v, j) => v - q[j]),
            )) /
          dot(m2, m2);
        const back = og.map((v, j) => v - f2 * m2[j]);
        expect(dot(g.n, back) - g.h).toBeLessThan(-1e-9);
        // the rotation of the plane (n_c, n_g) taking n_c to n_g agrees with the reflection on the directions tangent to H_c
        const c0 = dot(cl.n, g.n);
        const s0 = Math.sqrt(1 - c0 * c0);
        const e2 = Array.from(g.n, (v, j) => (v - c0 * cl.n[j]) / s0);
        const rot = (v: number[]) => {
          const p1 = dot(cl.n, v);
          const p2 = dot(e2, v);
          return v.map((x, j) => x + p1 * ((c0 - 1) * cl.n[j] + s0 * e2[j]) + p2 * (-s0 * cl.n[j] + (c0 - 1) * e2[j]));
        };
        const dir = Array.from(randomTangent(name, c, 1, rand));
        const viaRot = rot(dir);
        const viaRefl = dir.map((v, j) => v - ((2 * dot(m, dir)) / mm) * m[j]);
        viaRot.forEach((v, j) => expect(v).toBeCloseTo(viaRefl[j], 10));
      }
    }
  });
});

describe("the half-chain engine", () => {
  it("finds, by joins of half-chains, a path at least as short as every random geodesic walk, and its path is real (walking its direction ends at the other centre)", () => {
    const rand = mulberry(777);
    for (const name of ND_NAMES) {
      const C = getND(name);
      let checked = 0;
      for (let t = 0; t < 12; t++) {
        const s = ndCentreEdge(name) * (1.2 + 0.4 * rand());
        const c = C.pickCell(rand);
        const u = randomPoint(name, c, s, rand, 0.3);
        const x = C.toAmbient(c, u, s);
        const L = 0.3 + 0.69 * rand();
        const w = walkND(C, s, c, x, randomTangent(name, c, L, rand));
        if (w.hit) continue;
        checked++;
        const b = newNDBudget(5e6);
        const A = ndDevelop(C, s, c, x, 0.5 * L + 1e-9, b);
        const B = ndDevelop(C, s, w.cell, w.x, 0.5 * L + 1e-9, b);
        const j = ndJoinDistance(C, s, A, B, L + 1e-7, b, false);
        expect(j.dist, `${name} L = ${L}`).toBeLessThanOrEqual(L + 1e-9);
        const wv = new Float64Array(C.D);
        for (let r = 0; r < C.D; r++) wv[r] = B[j.b].a[r] - A[j.a].a[r];
        const z = ndDirectionAtAnchor(C, A, j.a, wv);
        const w2 = walkND(C, s, c, x, z);
        expect(w2.hit).toBe(false);
        expect(w2.cell).toBe(w.cell);
        expect(Math.hypot(...Array.from(w2.x, (v, r) => v - w.x[r]))).toBeLessThan(1e-8);
      }
      expect(checked, name).toBeGreaterThan(5);
    }
  }, 300_000);

  it("finds the hinge distance a random geodesic confirms: walking its nearest point just past it meets a hinge, just short of it does not", () => {
    const rand = mulberry(31);
    for (const name of ND_NAMES.filter((n) => /^(5|6|9)/.test(n))) {
      const C = getND(name);
      let checked = 0;
      for (let t = 0; t < 30 && checked < 3; t++) {
        const s = (name.includes("demicube") ? 1.0 : ndCentreEdge(name)) * (0.8 + 0.8 * rand());
        const c = C.pickCell(rand);
        const u = randomPoint(name, c, s, rand, 0.4);
        const x = C.toAmbient(c, u, s);
        const b = newNDBudget(2e4);
        try {
          const A = ndDevelop(C, s, c, x, 0.5, b);
          const h = ndHingeDistance(C, s, A, 0.5, b);
          if (!(h.dist < 0.5) || !h.x) continue;
          checked++;
          const wv = new Float64Array(C.D);
          for (let r = 0; r < C.D; r++) wv[r] = h.x[r] - A[h.node].a[r];
          const z = ndDirectionAtAnchor(C, A, h.node, wv);
          const zl = Math.hypot(...z);
          expect(
            walkND(
              C,
              s,
              c,
              x,
              z.map((v) => (v * h.dist * (1 + 1e-6)) / zl),
              null,
              1e-9,
            ).hit,
          ).toBe(true);
          expect(
            walkND(
              C,
              s,
              c,
              x,
              z.map((v) => (v * h.dist * (1 - 1e-6)) / zl),
              null,
              1e-9,
            ).hit,
          ).toBe(false);
        } catch (e) {
          if (!(e instanceof NDLimitError)) throw e;
        }
      }
      expect(checked, name).toBeGreaterThan(0);
    }
  }, 300_000);

  it("is cheap on packings: a few dozen chains per ball at reach 1/2", () => {
    for (const name of ND_NAMES) {
      const C = getND(name);
      const s = ndCentreEdge(name) * 1.3;
      const x = C.toAmbient(0, zeros(C.n), s);
      const b = newNDBudget(1e6);
      const nodes = ndDevelop(C, s, 0, x, 0.5, b);
      expect(nodes.length, name).toBeLessThan(400);
    }
  });
});

describe("known distances", () => {
  it("a cell's centre is at 0.9 x 1/2 from the hinges at 0.9 times the centre edge: s / sqrt(d (d-2)) on the simplex and orthoplex, s/2 on a demicube facet, s / sqrt 6 on the 5-demicube's", () => {
    for (const name of ND_NAMES) {
      const p = ndParse(name)!;
      const s = ndCentreEdge(name) * 0.9;
      const m = measureND(problem(name, 1), { s, pieces: [{ cell: 0, u: zeros(p.d - 1) }] });
      expect(m.clearance[0], name).toBeCloseTo(0.45, 12);
      // and above the centre edge nothing is reported: the clearance is more than the reach 1/2
      expect(
        measureND(problem(name, 1), { s: ndCentreEdge(name) * 1.05, pieces: [{ cell: 0, u: zeros(p.d - 1) }] })
          .clearance[0],
        name,
      ).toBe(Infinity);
    }
  });

  it("two simplex cells' centres are 2 r = s sqrt(2 / (n (n+1))) apart across their ridge, two adjacent demicube facets' centres s / sqrt 2", () => {
    for (const [d, s] of [
      [5, 3],
      [7, 4],
      [10, 6],
    ] as const) {
      const n = d - 1;
      const name = `${d}simplex` as NDName;
      const m = measureND(problem(name, 2), {
        s,
        pieces: [
          { cell: 0, u: zeros(n) },
          { cell: 1, u: zeros(n) },
        ],
      });
      expect(m.pairs).toHaveLength(1);
      expect(m.pairs[0].dist).toBeCloseTo(s * Math.sqrt(2 / (n * (n + 1))), 12);
    }
    for (const d of [6, 8]) {
      const name = `${d}demicube` as NDName;
      const s = 1.2;
      // cells 0 (x1 = -1/2) and 2 (x2 = -1/2) are adjacent at 90 degrees: the straight path across the ridge
      const m = measureND(problem(name, 2), {
        s,
        pieces: [
          { cell: 0, u: zeros(d - 1) },
          { cell: 2, u: zeros(d - 1) },
        ],
      });
      expect(m.pairs).toHaveLength(1);
      expect(m.pairs[0].dist).toBeCloseTo(s / Math.SQRT2, 12);
    }
  });

  it("matches the derived constructions: the simplex's d + 1 centres from sqrt(d (d-1)/2), the demicube's 2d from sqrt 2, two orthoplex neighbours from the simplex edge", () => {
    for (const name of ND_NAMES) {
      const p = ndParse(name)!;
      const rand = mulberry(5);
      const lay = ndExactLayout(name, p.family === "simplex" ? p.d + 1 : p.family === "orthoplex" ? 2 : 2 * p.d, rand)!;
      expect(lay).not.toBeNull();
      const cfg = {
        s: lay.s,
        pieces: Array.from(lay.st.cells, (cell, i) => ({
          cell,
          u: Array.from(lay.st.U.subarray(i * (p.d - 1), (i + 1) * (p.d - 1))),
        })),
      };
      expect(validateND(problem(name, cfg.pieces.length), cfg, { tol: 0 }).ok, name).toBe(true);
      // below the derived edge the pair of nearest centres overlaps
      const lo = { s: ndPairEdge(name) * (1 - 1e-6), pieces: cfg.pieces };
      expect(validateND(problem(name, cfg.pieces.length), lo, { tol: 0 }).ok, `${name} below`).toBe(false);
    }
  });
});

describe("the bounds", () => {
  it("are the larger of the volume bound and the one-ball bound, increasing in P", () => {
    for (const name of ND_NAMES) {
      const p = ndParse(name)!;
      let prev = 0;
      for (let P = 1; P <= ND_MAX_P; P++) {
        const b = ndBound(p.family, p.d, P);
        expect(b, `${name} ${P}`).toBeGreaterThanOrEqual(prev);
        prev = b;
        const one =
          p.family === "demicube"
            ? Math.sqrt((2 * (p.d - 2)) / ((p.d - 1) ** 2 - 5 * (p.d - 1) + 8))
            : Math.sqrt(p.d * (p.d - 2)) / 2;
        expect(b).toBeGreaterThanOrEqual(one * (1 - 1e-15));
        const vol = Math.pow((P * ndBallVolume(p.d - 1)) / ndBoundaryVolume(p.family, p.d), 1 / (p.d - 1));
        expect(b).toBeGreaterThanOrEqual(vol * (1 - 1e-12));
        expect(b).toBeLessThanOrEqual(Math.max(one, vol) * (1 + 1e-12));
        expect(ndBoundKind(p.family, p.d, P) === "volume").toBe(vol > one);
      }
    }
  });

  it("have the boundary volumes of the literature (unit edge)", () => {
    // 5-simplex: 6 facets, each a regular 4-simplex of edge 1, volume sqrt 5 / 96
    expect(ndBoundaryVolume("simplex", 5)).toBeCloseTo((6 * Math.sqrt(5)) / 96, 12);
    // 5-orthoplex: 32 facets of the same kind
    expect(ndBoundaryVolume("orthoplex", 5)).toBeCloseTo((32 * Math.sqrt(5)) / 96, 12);
    // 4D on 5D demicube: 10 16-cells of edge 1 (volume 1/3 each: the 4-cube of edge 1/sqrt 2 is 1/4, less 8 corners of 1/384... = 1/6)
    const sixteen = (1 - 8 / 24) / 4;
    expect(sixteen).toBeCloseTo(1 / 6, 14);
    expect(ndBoundaryVolume("demicube", 5)).toBeCloseTo(10 * sixteen + 16 * (Math.sqrt(5) / 96), 12);
  });

  it("are stated by the domain with their kind", () => {
    const e = boundaryDomain.bound({ fields: { polytope: "7demicube", piece: "ball" }, P: 3 });
    expect(e?.value).toBe(ndBound("demicube", 7, 3));
    expect(e?.label).toMatch(/one ball|volume/);
    expect(boundaryDomain.bound({ fields: { polytope: "7demicube", piece: "ball" }, P: 25 })).toBeUndefined();
    expect(boundaryDomain.bound({ fields: { polytope: "7simplex", piece: "cube" }, P: 3 })).toBeUndefined();
  });
});

describe("the float validator", () => {
  it("accepts the certifier's test packings at tolerance 0 and the server's, and refuses them shrunk by 1e-6", () => {
    expect(packings).toHaveLength(72);
    for (const r of packings) {
      const prob = problem(r.polytope, r.P);
      expect(validateND(prob, r.config, { tol: 0 }).ok, r.id).toBe(true);
      expect(
        ndServerRule(prob, r.config, ndBound(ndParse(r.polytope)!.family, ndParse(r.polytope)!.d, r.P)).valid,
        r.id,
      ).toBe(true);
      const f = 1 - 1e-6;
      const shrunk = {
        s: r.config.s * f,
        pieces: r.config.pieces.map((q) => ({ cell: q.cell, u: q.u.map((x) => x * f) })),
      };
      expect(validateND(prob, shrunk, { tol: 1e-9 }).ok, `${r.id} shrunk`).toBe(false);
    }
  });

  it("refuses malformed input with reasons and never throws", () => {
    const prob = problem("7orthoplex", 2);
    const ok = {
      s: 5,
      pieces: [
        { cell: 0, u: zeros(6) },
        { cell: 127, u: zeros(6) },
      ],
    };
    expect(validateND(prob, ok).ok).toBe(true);
    const bad: unknown[] = [
      null,
      [],
      {},
      { s: "5", pieces: ok.pieces },
      { s: NaN, pieces: ok.pieces },
      { s: -1, pieces: ok.pieces },
      { s: 5, pieces: [ok.pieces[0]] },
      { s: 5, pieces: [ok.pieces[0], { cell: 128, u: zeros(6) }] },
      { s: 5, pieces: [ok.pieces[0], { cell: -1, u: zeros(6) }] },
      { s: 5, pieces: [ok.pieces[0], { cell: 1.5, u: zeros(6) }] },
      { s: 5, pieces: [ok.pieces[0], { cell: 1, u: zeros(5) }] },
      { s: 5, pieces: [ok.pieces[0], { cell: 1, u: [0, 0, 0, 0, 0, Infinity] }] },
      { s: 5, pieces: [ok.pieces[0], { cell: 1, u: zeros(6), q: [1, 0, 0, 0] }] },
      { s: 5, pieces: [ok.pieces[0], { cell: 1, u: [0, 0, 0, 0, 0, 99] }] },
    ];
    for (const b of bad) {
      const v = validateND(prob, b);
      expect(v.ok, JSON.stringify(b)).toBe(false);
      expect(v.reasons.length).toBeGreaterThan(0);
    }
    // a sparse array: a hole is rejected, not read as 0
    const holes = [0, 0, , 0, 0, 0] as unknown as number[];
    expect(validateND(prob, { s: 5, pieces: [ok.pieces[0], { cell: 1, u: holes }] }).ok).toBe(false);
    expect(validateND({ polytope: "7orthoplex", piece: "cube", P: 2 } as never, ok).ok).toBe(false);
    expect(validateND({ polytope: "8cell", piece: "ball", P: 2 } as never, ok).ok).toBe(false);
  });

  it("fails closed at the step budget: a check that runs out is refused with a reason, never accepted", () => {
    const r = packings.find((q) => q.polytope === "6demicube" && q.P === 4)!;
    const prob = problem(r.polytope, r.P);
    const full = { steps: 0 };
    expect(validateND(prob, r.config, { tol: 1e-9, stats: full }).ok).toBe(true);
    expect(full.steps).toBeGreaterThan(100);
    for (const limit of [1, 5, 50, Math.floor(full.steps / 2), full.steps - 1]) {
      const v = validateND(prob, r.config, { tol: 1e-9, maxWork: limit });
      expect(v.ok, `limit ${limit}`).toBe(false);
      expect(v.reasons[0]).toMatch(/work limit/);
    }
    expect(validateND(prob, r.config, { tol: 1e-9, maxWork: full.steps }).ok).toBe(true);
  });

  it("refuses a centre outside its cell, a ball over a hinge, and balls too close across a ridge, and accepts them just clear", () => {
    const name: NDName = "6demicube";
    const C = getND(name);
    const s = 1.4;
    // outside its cell: a demicube facet's coordinates are |u_k| <= s/2 ... the point beyond a ridge
    const out = { s, pieces: [{ cell: 0, u: [s, 0, 0, 0, 0] }] };
    expect(validateND(problem(name, 1), out).reasons[0]).toMatch(/outside its anchor cell/);
    // a ball at the centre of a facet, then one moved toward a hinge until its clearance is 0.5 + 1e-3 (valid) and 0.5 - 1e-3 (invalid)
    const clearance = (t: number) =>
      measureND(problem(name, 1), { s, pieces: [{ cell: 0, u: [t, t, 0, 0, 0] }] }).clearance[0];
    let lo = 0;
    let hi = s / 2;
    for (let i = 0; i < 60; i++) {
      const mid = (lo + hi) / 2;
      if (clearance(mid) > 0.5) lo = mid;
      else hi = mid;
    }
    expect(clearance(lo)).toBeGreaterThan(0.5 - 1e-9);
    const inside = { s, pieces: [{ cell: 0, u: [lo * 0.999, lo * 0.999, 0, 0, 0] }] };
    const over = { s, pieces: [{ cell: 0, u: [hi * 1.002, hi * 1.002, 0, 0, 0] }] };
    expect(validateND(problem(name, 1), inside, { tol: 0 }).ok).toBe(true);
    expect(validateND(problem(name, 1), over, { tol: 0 }).ok).toBe(false);
    // two balls on either side of a ridge centre, 1 + 1e-6 and 1 - 1e-6 apart along the common perpendicular
    const cl = C.cell(0);
    void cl;
    const near = (gap: number) => {
      // cells 0 (x1 = -1/2) and 2 (x2 = -1/2) are adjacent at 90 degrees; the facets' planes are at s / (2 sqrt 2) from the centre of the
      // container (unit edge), and the straight path across the ridge between points at a from it in cell 0 and in cell 2 has length 2a:
      // a = (1 + gap) / 2 (each ball then touches the hinge, at distance a)
      const a = (1 + gap) / 2;
      const h = s / (2 * Math.SQRT2);
      return {
        s,
        pieces: [
          // in cell 0 (x1 = -h) the local coordinates are x2..x6: the point at x2 = -h + a
          { cell: 0, u: [-h + a, 0, 0, 0, 0] },
          // in cell 2 (x2 = -h) the local coordinates are x1, x3..x6: x1 = -h + a
          { cell: 2, u: [-h + a, 0, 0, 0, 0] },
        ],
      };
    };
    expect(validateND(problem(name, 2), near(1e-6), { tol: 0 }).ok).toBe(true);
    expect(validateND(problem(name, 2), near(-1e-6), { tol: 1e-9 }).ok).toBe(false);
    expect(validateND(problem(name, 2), near(-1e-6), { tol: 1e-9 }).reasons[0]).toMatch(/overlap/);
  });
});

describe("the server rule, the registry and the heavy queue", () => {
  it("parses the ids and lists 18 containers x 24 problems", () => {
    expect(ND_NAMES).toHaveLength(18);
    for (const name of ND_NAMES) {
      expect(isContainerName(name)).toBe(true);
      expect(CONTAINERS).toContain(name);
      expect(CONTAINER_MAX_P[name]).toBe(24);
      expect(parseBoundaryId(`boundary/${name}/ball/24`)).toEqual({ polytope: name, piece: "ball", P: 24 });
      expect(typeof parseBoundaryId(`boundary/${name}/cube/3`)).toBe("string");
      expect(typeof parseBoundaryId(`boundary/${name}/ball/07`)).toBe("string");
    }
    expect(boundaryDomain.info.problems).toBeGreaterThan(432);
    const n = boundaryDomain.info.fields.find((f) => f.key === "polytope")!;
    for (const name of ND_NAMES) expect(n.values).toContain(name);
    expect(
      boundaryDomain.sample(100_000).filter((id) => /\/(\d+)(simplex|orthoplex|demicube)\//.test(id)),
    ).toHaveLength(432);
  });

  it("refuses P above 24, wrong cells, a rotation and a size below the bound, and accepts the fixtures", () => {
    const name: NDName = "8simplex";
    const p = ndParse(name)!;
    const b = ndBound(p.family, p.d, 2);
    const okc = {
      s: 5.3,
      pieces: [
        { cell: 0, u: zeros(7) },
        { cell: 8, u: zeros(7) },
      ],
    };
    expect(ndServerRule(problem(name, 2), okc, b).valid).toBe(true);
    expect(ndServerRule(problem(name, 25), okc, b).valid).toBe(false);
    expect(ndServerRule(problem(name, 2), { ...okc, pieces: [okc.pieces[0], { cell: 9, u: zeros(7) }] }, b).valid).toBe(
      false,
    );
    expect(
      ndServerRule(problem(name, 2), { ...okc, pieces: [okc.pieces[0], { cell: 8, u: zeros(7), q: [1, 0, 0, 0] }] }, b)
        .valid,
    ).toBe(false);
    expect(ndServerRule(problem(name, 2), { ...okc, s: b * (1 - 1e-6) }, b).valid).toBe(false);
  });

  it("puts the checks that may take 50 ms or more in the heavy queue, from the measured per-ball costs", () => {
    expect(boundaryIsHeavy({ polytope: "8demicube", piece: "ball", P: 1 })).toBe(true);
    expect(boundaryIsHeavy({ polytope: "7simplex", piece: "ball", P: 24 })).toBe(false);
    expect(boundaryIsHeavy({ polytope: "10orthoplex", piece: "ball", P: 24 })).toBe(false);
    expect(boundaryIsHeavy({ polytope: "10demicube", piece: "ball", P: 12 })).toBe(false);
    expect(boundaryIsHeavy({ polytope: "10demicube", piece: "ball", P: 13 })).toBe(true);
    expect(BOUNDARY_HEAVY_MS).toBe(50);
  });
});

describe("the kernel and the search", () => {
  const kernel = (name: NDName, P: number) => boundaryKernel({ polytope: name, piece: "ball", P });
  it("passes the kernel contract on a container of each family and dimension range", () => {
    for (const [name, P] of [
      ["5simplex", 3],
      ["7simplex", 6],
      ["10simplex", 9],
      ["6orthoplex", 4],
      ["9orthoplex", 6],
      ["5demicube", 5],
      ["7demicube", 3],
      ["10demicube", 4],
    ] as const)
      expect(
        checkKernel(() => kernel(name, P), { seed: P, states: 2 }),
        `${name} ${P}`,
      ).toEqual([]);
  }, 600_000);

  it("seeds valid layouts, exact where one is derived and sequential otherwise, each with no energy", () => {
    const rng = mulberry(4);
    for (const [name, P] of [
      ["6simplex", 3],
      ["6simplex", 20],
      ["8orthoplex", 2],
      ["8orthoplex", 11],
      ["6demicube", 8],
      ["6demicube", 24],
    ] as const) {
      const k = kernel(name, P);
      const seeds = k.seeds!(rng, 2);
      expect(seeds).toHaveLength(2);
      for (const q of seeds) {
        expect(k.energy(q), `${name} ${P}`).toBeLessThan(1e-20);
        expect(k.validate(k.toConfig(q)).reasons).toEqual([]);
      }
    }
    const lay = ndSequentialLayout("7orthoplex", 10, mulberry(2));
    expect(
      validateND(
        problem("7orthoplex", 10),
        {
          s: lay.s,
          pieces: Array.from(lay.st.cells, (cell, i) => ({
            cell,
            u: Array.from(lay.st.U.subarray(i * 6, (i + 1) * 6)),
          })),
        },
        { tol: 0 },
      ).ok,
    ).toBe(true);
  }, 300_000);

  it("runs the boundary search: seeded, reproducible, valid at tolerance 0, never below the bound", () => {
    for (const [name, P] of [
      ["5simplex", 4],
      ["8orthoplex", 3],
      ["6demicube", 6],
    ] as const) {
      const a = searchBoundary({ polytope: name, piece: "ball", P }, { seed: 3, budget: { evals: 3000 } });
      const b = searchBoundary({ polytope: name, piece: "ball", P }, { seed: 3, budget: { evals: 3000 } });
      expect(a.config).not.toBeNull();
      expect(a.score).toBe(b.score);
      expect(validateBoundary({ polytope: name, piece: "ball", P }, a.config, { tol: 0 }).ok).toBe(true);
      const p = ndParse(name)!;
      expect(a.score!).toBeGreaterThanOrEqual(ndBound(p.family, p.d, P));
    }
  }, 300_000);
});
