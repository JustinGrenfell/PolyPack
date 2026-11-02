/**
 * The engine's distances on configurations of balls on the simplex, orthoplex and demicube (d = 5..10), for python/boundary/test_ndpoly.py
 * to compare with the exact certifier (docs/spaces/BOUNDARY.md 4D.6). Seeded random layouts of 1 to 3 balls on each of the 18 containers:
 * the first ball at a random point of a random cell (nearer the cell's boundary on every third case), each next ball at the end of a random
 * geodesic of length 0.7 to 1.2 from a random earlier ball (so that most pairs cross ridges and many are close to 1 apart). Per case: the
 * engine's clearance of each centre, the pair distances below 1 and the loops below 1 (measureND), and its verdict at tolerance 0.
 * Written to python/boundary/fixtures/nd_crosscheck.json.
 *
 *   npx tsx packages/boundary/scripts/gen-nd-crosscheck.ts
 */
import { writeFileSync } from "node:fs";
import { fileURLToPath } from "node:url";
import { ND_NAMES, type NDName, getND, measureND, ndCentreEdge, ndParse, validateND, walkND } from "../src/index";
import { mulberry32 } from "../src/search/random";

const out = fileURLToPath(new URL("../../../python/boundary/fixtures/nd_crosscheck.json", import.meta.url));
const PER = Number(process.argv[2] ?? 4);
const rand = mulberry32(20261102);
const gauss = () => Math.sqrt(-2 * Math.log(1 - rand())) * Math.cos(2 * Math.PI * rand());

interface Case {
  polytope: string;
  s: number;
  pieces: { cell: number; u: number[] }[];
  clearance: (number | null)[];
  pairs: { i: number; j: number; dist: number | null }[];
  loops: (number | null)[];
  valid0: boolean;
  steps: number;
}
const cases: Case[] = [];
const fin = (x: number) => (Number.isFinite(x) ? x : null);

for (const name of ND_NAMES) {
  const C = getND(name);
  const p = ndParse(name)!;
  for (let t = 0; t < PER; t++) {
    const s = ndCentreEdge(name) * (1.0 + 0.7 * rand());
    const P = 1 + Math.floor(rand() * 3);
    const cell = C.pickCell(rand);
    const cl = C.cell(cell);
    const spread = (t % 3 === 2 ? 0.5 : 0.2) / Math.sqrt(C.n);
    let u = Array.from({ length: C.n }, () => gauss() * spread * s);
    let x = C.toAmbient(cell, u, s);
    const inside = (xx: Float64Array) => {
      for (let k = 0; k < cl.nbr.length; k++) {
        let v = -s * cl.hn[k];
        for (let j = 0; j < C.D; j++) v += cl.nn[k * C.D + j] * xx[j];
        if (v > 0) return false;
      }
      return true;
    };
    let lam = 1;
    const u0 = u;
    while (!inside(x) && lam > 1e-3) {
      lam *= 0.8;
      u = u0.map((v) => v * lam);
      x = C.toAmbient(cell, u, s);
    }
    if (!inside(x)) continue;
    const pieces = [{ cell, u }];
    const xs = [{ cell, x }];
    for (let q = 1; q < P; q++) {
      const from = xs[Math.floor(rand() * xs.length)];
      const fc = C.cell(from.cell);
      for (let tries = 0; tries < 20; tries++) {
        const dv = new Float64Array(C.D);
        for (let k = 0; k < C.n; k++) {
          const g = gauss();
          for (let j = 0; j < C.D; j++) dv[j] += g * fc.frame[k][j];
        }
        let nn = 0;
        for (let j = 0; j < C.D; j++) nn += dv[j] * dv[j];
        const L = 0.7 + 0.5 * rand();
        for (let j = 0; j < C.D; j++) dv[j] *= L / Math.sqrt(nn);
        const w = walkND(C, s, from.cell, from.x, dv);
        if (w.hit) continue;
        pieces.push({ cell: w.cell, u: Array.from(C.toLocal(w.cell, w.x, s)) });
        xs.push({ cell: w.cell, x: w.x });
        break;
      }
    }
    const config = { s, pieces };
    const problem = { polytope: name, piece: "ball" as const, P: pieces.length };
    try {
      const m = measureND(problem, config);
      const v = validateND(problem, config, { tol: 0 });
      cases.push({
        polytope: name,
        s,
        pieces,
        clearance: m.clearance.map(fin),
        pairs: m.pairs.map((r) => ({ i: r.i, j: r.j, dist: fin(r.dist) })),
        loops: m.loops.map(fin),
        valid0: v.ok,
        steps: m.steps,
      });
    } catch (e) {
      console.error("skipped", name, String(e).slice(0, 80));
    }
  }
  void p;
}
// ---- adversarial cases: pairs near 1 apart, balls near 1/2 from a hinge, balls on a ridge
const adv: (Case & { kind: string; delta: number })[] = [];
const record = (name: NDName, s: number, pieces: Case["pieces"], kind: string, delta: number) => {
  const problem = { polytope: name, piece: "ball" as const, P: pieces.length };
  const config = { s, pieces };
  try {
    const m = measureND(problem, config);
    const v = validateND(problem, config, { tol: 0 });
    adv.push({
      kind,
      delta,
      polytope: name,
      s,
      pieces,
      clearance: m.clearance.map(fin),
      pairs: m.pairs.map((r) => ({ i: r.i, j: r.j, dist: fin(r.dist) })),
      loops: m.loops.map(fin),
      valid0: v.ok,
      steps: m.steps,
    });
  } catch (e) {
    console.error("adv skipped", name, kind, String(e).slice(0, 80));
  }
};
const randomDir = (C: ReturnType<typeof getND>, cell: number, len: number): Float64Array => {
  const fc = C.cell(cell);
  const dv = new Float64Array(C.D);
  for (let k = 0; k < C.n; k++) {
    const g = gauss();
    for (let j = 0; j < C.D; j++) dv[j] += g * fc.frame[k][j];
  }
  let nn = 0;
  for (let j = 0; j < C.D; j++) nn += dv[j] * dv[j];
  for (let j = 0; j < C.D; j++) dv[j] *= len / Math.sqrt(nn);
  return dv;
};
for (const name of ND_NAMES) {
  const C = getND(name);
  const s = ndCentreEdge(name) * 1.25;
  // pairs: ball 2 at the end of a geodesic of length 1 + delta from a ball near its cell's centre
  for (const delta of [1e-3, 1e-6, 1e-9, 1e-12, -1e-12, -1e-9, -1e-6]) {
    for (let tries = 0; tries < 40; tries++) {
      const cell = C.pickCell(rand);
      const cl0 = C.cell(cell);
      const u0 = Array.from({ length: C.n }, () => (gauss() * 0.1 * s) / Math.sqrt(C.n));
      // shrunk toward the cell's centre until the centre is inside the cell
      let u = u0;
      let x = C.toAmbient(cell, u, s);
      for (let lam = 0.8, k = 0; k < 40; k++, lam *= 0.8) {
        let ok = true;
        for (let q = 0; q < cl0.nbr.length && ok; q++) {
          let v = -s * cl0.hn[q];
          for (let j = 0; j < C.D; j++) v += cl0.nn[q * C.D + j] * x[j];
          if (v > 0) ok = false;
        }
        if (ok) break;
        u = u0.map((t) => t * lam);
        x = C.toAmbient(cell, u, s);
      }
      const w = walkND(C, s, cell, x, randomDir(C, cell, 1 + delta));
      if (w.hit || w.cells.length < 2) continue;
      record(
        name,
        s,
        [
          { cell, u },
          { cell: w.cell, u: Array.from(C.toLocal(w.cell, w.x, s)) },
        ],
        "pair",
        delta,
      );
      break;
    }
  }
  // hinge: the centre of a cell moved away from its nearest hinge until the clearance is 1/2 + delta (bisection on the engine's own measure)
  for (const delta of [1e-3, 1e-9, 1e-12, -1e-12, -1e-9, -1e-3]) {
    const cell = name.endsWith("demicube") ? 0 : C.pickCell(rand);
    const dirU = Array.from({ length: C.n }, () => gauss());
    const nm = Math.sqrt(dirU.reduce((a, v) => a + v * v, 0));
    const at = (t: number) => dirU.map((v) => (v / nm) * t);
    const clear = (t: number) =>
      measureND({ polytope: name, piece: "ball", P: 1 }, { s, pieces: [{ cell, u: at(t) }] }).clearance[0];
    // clearance falls from its centre value as the ball moves out: find t in [0, tmax] with clearance 1/2 + delta
    let lo = 0;
    let hi = 0.6 * s;
    if (!(clear(lo) > 0.5 + delta)) continue;
    // keep the point inside the cell: shrink hi until it is
    const cl = C.cell(cell);
    const inCell = (t: number) => {
      const xx = C.toAmbient(cell, at(t), s);
      for (let k = 0; k < cl.nbr.length; k++) {
        let v = -s * cl.hn[k];
        for (let j = 0; j < C.D; j++) v += cl.nn[k * C.D + j] * xx[j];
        if (v > 0) return false;
      }
      return true;
    };
    while (!inCell(hi)) hi *= 0.8;
    if (!(clear(hi) < 0.5 + delta)) continue;
    for (let it = 0; it < 60; it++) {
      const mid = 0.5 * (lo + hi);
      if (clear(mid) > 0.5 + delta) lo = mid;
      else hi = mid;
    }
    record(name, s, [{ cell, u: at(lo) }], "hinge", delta);
  }
  // a ball on a ridge (exactly, by projecting onto the plane of a neighbour) and a second one 0.9 and 1.05 away across it
  for (const dist of [0.9, 1.05]) {
    const cell = C.pickCell(rand);
    const cl = C.cell(cell);
    const k = Math.floor(rand() * cl.nbr.length);
    // the ridge point of neighbour k, nudged toward the cell's centre
    const q = Float64Array.from(cl.rp.subarray(k * C.D, k * C.D + C.D), (v) => v * s);
    const o = Float64Array.from(cl.o, (v) => v * s);
    const x = Float64Array.from(q, (v, j) => v + 0.3 * (o[j] - v));
    // onto the plane of neighbour k exactly
    let e = -s * cl.hn[k];
    for (let j = 0; j < C.D; j++) e += cl.nn[k * C.D + j] * x[j];
    const c1 = cl.nn.subarray(k * C.D, k * C.D + C.D).reduce((a, v, j) => a + v * cl.n[j], 0);
    const nu = Float64Array.from(cl.nn.subarray(k * C.D, k * C.D + C.D), (v, j) => v - c1 * cl.n[j]);
    for (let j = 0; j < C.D; j++) x[j] -= (e / cl.nu2[k]) * nu[j];
    const u = Array.from(C.toLocal(cell, x, s));
    const w = walkND(C, s, cell, x, randomDir(C, cell, dist));
    if (w.hit) continue;
    record(
      name,
      s,
      [
        { cell, u },
        { cell: w.cell, u: Array.from(C.toLocal(w.cell, w.x, s)) },
      ],
      "ridge",
      dist,
    );
  }
}
writeFileSync(out, JSON.stringify({ about: "gen-nd-crosscheck.ts", cases, adversarial: adv }, null, 0) + "\n");
console.log(`wrote ${cases.length} cases and ${adv.length} adversarial cases`);
