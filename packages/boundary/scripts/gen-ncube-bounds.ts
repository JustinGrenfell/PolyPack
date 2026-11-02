/**
 * The engine's lower bounds on s for the d-cube problems (a list of P per d and piece), for python/boundary/test_ncube.py to check
 * exactly. Written to python/boundary/fixtures/ncube_bounds.json.
 *
 *   npx tsx packages/boundary/scripts/gen-ncube-bounds.ts
 */
import { writeFileSync } from "node:fs";
import { fileURLToPath } from "node:url";
import { NCUBE_MAX_P, ncubeBound } from "../src/index";

const out = fileURLToPath(new URL("../../../python/boundary/fixtures/ncube_bounds.json", import.meta.url));
const bounds: { d: number; piece: string; P: number; value: number }[] = [];
for (let d = 5; d <= 10; d++) {
  const max = NCUBE_MAX_P[`${d}cube` as keyof typeof NCUBE_MAX_P];
  const Ps = new Set<number>(
    [1, 2, 3, d, 2 * d, 2 * d + 1, 3 * d, 4 * d, 10 * d, 2 * d * 2 ** (d - 1), max, max - 1, 50, 100, 160, 384].filter(
      (P) => P >= 1 && P <= max,
    ),
  );
  for (const P of [...Ps].sort((a, b) => a - b))
    for (const piece of ["ball", "cube"] as const) bounds.push({ d, piece, P, value: ncubeBound(d, piece, P) });
}
writeFileSync(
  out,
  JSON.stringify({
    note: "written by packages/boundary/scripts/gen-ncube-bounds.ts; read by python/boundary/test_ncube.py",
    bounds,
  }),
);
console.log("bounds", bounds.length);
