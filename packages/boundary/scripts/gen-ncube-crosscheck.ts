/**
 * The engine's verdicts on d-cube configurations (d = 5..10), for python/boundary/test_ncube.py to compare with the exact
 * certifier: the small test packings of fixtures/ncube_packings.json scaled by 1.2, 1, 0.9999 and 0.97 (valid, touching or
 * barely valid, and overlapping variants), and seeded random layouts of 2 to 5 pieces, half of them put near the ridges.
 * Written to python/boundary/fixtures/ncube_crosscheck.json.
 *
 *   npx tsx packages/boundary/scripts/gen-ncube-crosscheck.ts
 */
import { readFileSync, writeFileSync } from "node:fs";
import { fileURLToPath } from "node:url";
import { NCUBE_NAMES, validateNCube } from "../src/index";
import { mulberry32 } from "../src/search/random";

const fixtures = fileURLToPath(new URL("../../../python/boundary/fixtures/", import.meta.url));
const packings: {
  id: string;
  piece: "ball" | "cube";
  P: number;
  polytope: string;
  config: { s: number; pieces: { cell: number; u: number[] }[] };
}[] = JSON.parse(readFileSync(fixtures + "ncube_packings.json", "utf8")).packings;

interface Case {
  d: number;
  piece: "ball" | "cube";
  s: number;
  pieces: { cell: number; u: number[] }[];
  tol9: boolean;
  tol0: boolean;
  reason: string | null;
  source: string;
}
const cases: Case[] = [];
const add = (d: number, piece: "ball" | "cube", s: number, pieces: Case["pieces"], source: string) => {
  const problem = { polytope: `${d}cube` as (typeof NCUBE_NAMES)[number], piece, P: pieces.length };
  const a = validateNCube(problem, { s, pieces }, { tol: 1e-9 });
  const b = validateNCube(problem, { s, pieces }, { tol: 0 });
  cases.push({ d, piece, s, pieces, tol9: a.ok, tol0: b.ok, reason: b.reasons[0] ?? null, source });
};

for (const p of packings)
  for (const f of [1.2, 1, 0.9999, 0.97]) {
    const d = Number(p.polytope.slice(0, -4));
    add(
      d,
      p.piece,
      p.config.s * f,
      p.config.pieces.map((q) => ({ cell: q.cell, u: q.u.map((x) => x * f) })),
      `${p.id} x ${f}`,
    );
  }

const rand = mulberry32(2024);
for (let d = 5; d <= 10; d++) {
  const n = d - 1;
  for (const piece of ["ball", "cube"] as const)
    for (let t = 0; t < 40; t++) {
      const P = 2 + Math.floor(rand() * 4);
      const s = piece === "ball" ? 0.8 + rand() * 1.4 : 1 + rand() * 1.2;
      const h = s / 2;
      const pieces: Case["pieces"] = [];
      for (let i = 0; i < P; i++) {
        // near the ridges: some coordinates within 0.45 of a face; each piece valid alone (a clear ball, a cube crossing at most one axis)
        let u: number[] = [];
        for (let a = 0; a < 500; a++) {
          u = Array.from({ length: n }, () =>
            rand() < 0.3 ? (rand() < 0.5 ? -1 : 1) * (h - rand() * 0.45) : (2 * rand() - 1) * h * 0.8,
          ).map((x) => Math.max(-h, Math.min(h, x)));
          const slack = u.map((x) => h - Math.abs(x)).sort((p, q) => p - q);
          if (
            piece === "ball"
              ? Math.hypot(slack[0], slack[1]) >= 0.5
              : u.filter((x) => Math.abs(x) + 0.5 > h).length <= 1
          )
            break;
        }
        pieces.push({ cell: Math.floor(rand() * 2 * d), u });
      }
      add(d, piece, s, pieces, `random ${d} ${piece} ${t}`);
    }
}
const valid = cases.filter((c) => c.tol0).length;
console.log("cases", cases.length, "valid at tol 0", valid, "valid at 1e-9", cases.filter((c) => c.tol9).length);
writeFileSync(
  fixtures + "ncube_crosscheck.json",
  JSON.stringify({
    note: "written by packages/boundary/scripts/gen-ncube-crosscheck.ts; read by python/boundary/test_ncube.py",
    cases,
  }),
);
