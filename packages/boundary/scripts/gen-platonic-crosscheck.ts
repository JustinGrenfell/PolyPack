/**
 * Cross-check fixture of the float validator against the exact certifier for Platonic solids (BOUNDARY.md section 4C.5): for each
 * container and solid, two solids at random cells, random turns, and (a third of the time) centres on a ridge, scaled about the centres of the
 * cells to the critical scale k* at which the engine's validator (tolerance 0) stops accepting them (found by bisection: validity is monotone in
 * the scale). At k* (1 + d) and k* (1 - d), d = 1e-5 and 1e-9, the engine's verdicts at tolerance 0 and 1e-9 are recorded; the Python test
 * (python/boundary/test_platonic.py) must reach the same verdict exactly, whatever contact it is: two solids touching, an edge of the
 * polytope touching a solid, a solid touching itself around a sharp edge, in one cell or across ridges.
 *
 *   npx tsx packages/boundary/scripts/gen-platonic-crosscheck.ts
 */
import { writeFileSync } from "node:fs";
import { fileURLToPath } from "node:url";
import type { BoundaryConfig, PolytopeName, PolytopePieceName } from "../src/index";
import { CONTAINER_MAX_P, cellExcess, getPolytope, matrixToQuaternion, pieceBounds, validateBoundary } from "../src/index";
import { randomCellPoint } from "../src/search/chart";
import { mulberry32, randInt } from "../src/search/random";
import { randomRotation } from "../src/search/rot";
import { SolidModel } from "../src/search/solids";

const out = fileURLToPath(new URL("../../../python/boundary/fixtures/platonic_crosscheck.json", import.meta.url));
const NAMES: PolytopeName[] = ["tesseract", "5cell", "16cell", "24cell", "120cell", "600cell"];
const PIECES: PolytopePieceName[] = ["tetrahedron", "cube", "octahedron", "dodecahedron", "icosahedron"];
const rand = mulberry32(2026);
const entries: Record<string, unknown>[] = [];

const scale = (c: BoundaryConfig, k: number): BoundaryConfig => ({
  s: c.s * k,
  pieces: c.pieces.map((p) => ({ cell: p.cell, u: p.u.map((x) => x * k), ...(p.q ? { q: p.q } : {}) })),
});

for (const name of NAMES)
  for (const piece of PIECES) {
    if (name === "tesseract" && piece === "cube") continue;
    const poly = getPolytope(name);
    const model = new SolidModel(poly, piece, 2);
    const problem = { polytope: name, piece, P: 2 };
    const check = (c: BoundaryConfig, tol: number) => validateBoundary(problem, c, { tol, maxChains: 4096 });
    let made = 0;
    let onRidge = false;
    for (let attempt = 0; attempt < 120 && made < 2; attempt++) {
      // the first entry: random points, every other attempt with a centre on a ridge; the second: the centres of two adjacent cells
      // (the solids meet across their common ridge, or around its edges)
      const s0 = model.startEdge() * 1.2;
      const pieces = [];
      onRidge = false;
      const cell0 = randInt(rand, poly.cells.length);
      for (let i = 0; i < 2; i++) {
        const cell = made === 1 ? (i === 0 ? cell0 : poly.cells[cell0].facets[randInt(rand, poly.cells[cell0].facets.length)].neighbor) : randInt(rand, poly.cells.length);
        let u = made === 1 ? [0, 0, 0] : randomCellPoint(poly, cell, s0, rand, 1);
        if (made === 0 && attempt % 2 === 0 && i === 0) {
          // onto a random facet plane of the cell (a ridge), kept inside the closed cell
          const f = poly.cells[cell].facets[randInt(rand, poly.cells[cell].facets.length)];
          const e = f.inside.n[0] * u[0] + f.inside.n[1] * u[1] + f.inside.n[2] * u[2] - f.inside.h * s0;
          const v = u.map((x, k) => x - e * f.inside.n[k]);
          if (cellExcess(poly, cell, v, s0) <= 0) {
            u = v;
            onRidge = true;
          }
        }
        pieces.push({ cell, u, q: matrixToQuaternion(randomRotation(rand)) });
      }
      const base: BoundaryConfig = { s: s0, pieces };
      // The tolerance of the comparison: 0, except when a centre lies on a ridge (its float coordinates can then be outside the cell by a rounding
      // unit, which tolerance 0 refuses whatever the contact): 1e-9, the server's.
      const tol = onRidge ? 1e-9 : 0;
      const valid = (c: BoundaryConfig, t: number) => check(c, t);
      // grow the scale until valid, then bisect down to the critical one
      let hi = 1;
      let tries = 0;
      while (!valid(scale(base, hi), tol).ok && tries++ < 12) hi *= 1.25;
      if (tries > 12 || hi > 3) continue;
      let lo = hi;
      let okLo = false;
      for (let k = 0; k < 40 && !okLo; k++) {
        lo *= 0.9;
        if (!valid(scale(base, lo), tol).ok) okLo = true;
      }
      if (!okLo) continue;
      for (let it = 0; it < 60; it++) {
        const mid = 0.5 * (lo + hi);
        if (valid(scale(base, mid), tol).ok) hi = mid;
        else lo = mid;
      }
      const kStar = hi;
      const verdicts: Record<string, unknown>[] = [];
      // the 1e-9 cases only where the comparison is at tolerance 0 (at 1e-9 the engine's own threshold is that depth)
      for (const d of entries.length % 4 === 0 && tol === 0 ? [1e-5, 1e-9] : [1e-5])
        for (const sign of [1, -1]) {
          const cfg = scale(base, kStar * (1 + sign * d));
          const v = valid(cfg, tol);
          verdicts.push({ d, sign, config: cfg, engine: v.ok, reason: v.reasons[0] ?? null });
        }
      // keep it only when the critical contact is real: the engine accepts the loose side and refuses the tight side at 1e-5
      const [plus, minus] = verdicts;
      if (!(plus.engine && !minus.engine)) continue;
      entries.push({ polytope: name, piece, P: 2, ridge: onRidge, adjacent_centres: made === 1, tol, cases: verdicts });
      made++;
      console.log(name, piece, "k*", kStar.toFixed(6), "s*", (s0 * kStar).toFixed(5), minus.reason ? String(minus.reason).slice(0, 70) : "");
    }
  }
// the engine's lower bounds of every Platonic problem, for the exact check in the Python test
const bounds: Record<string, number> = {};
for (const name of NAMES)
  for (const piece of PIECES) {
    if (name === "tesseract" && piece === "cube") continue;
    for (let P = 1; P <= CONTAINER_MAX_P[name]; P++) bounds[`boundary/${name}/${piece}/${P}`] = pieceBounds(name, piece, P).value;
  }
writeFileSync(out, JSON.stringify({ note: "generated by packages/boundary/scripts/gen-platonic-crosscheck.ts (seed 2026)", entries, bounds }));
console.log("wrote", entries.length, "entries");
