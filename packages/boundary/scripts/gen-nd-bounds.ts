/**
 * The engine's lower bounds on s for the simplex, orthoplex and demicube problems (every P = 1..24 on every container), for
 * python/boundary/test_ndpoly.py to check exactly. Written to python/boundary/fixtures/nd_bounds.json.
 *
 *   npx tsx packages/boundary/scripts/gen-nd-bounds.ts
 */
import { writeFileSync } from "node:fs";
import { fileURLToPath } from "node:url";
import { ND_MAX_P, ND_NAMES, ndBound, ndParse } from "../src/index";

const out = fileURLToPath(new URL("../../../python/boundary/fixtures/nd_bounds.json", import.meta.url));
const bounds: { polytope: string; P: number; value: number }[] = [];
for (const name of ND_NAMES) {
  const p = ndParse(name)!;
  for (let P = 1; P <= ND_MAX_P; P++) bounds.push({ polytope: name, P, value: ndBound(p.family, p.d, P) });
}
writeFileSync(
  out,
  JSON.stringify({
    note: "written by packages/boundary/scripts/gen-nd-bounds.ts; read by python/boundary/test_ndpoly.py",
    bounds,
  }),
);
console.log("bounds", bounds.length);
