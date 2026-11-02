/**
 * Writes python/boundary/fixtures/containers.json: the conventions of the six containers as the engine builds them (vertex
 * lists, and for every k-th cell its vertex ids, origin, frame and normal), which python/boundary/test_general.py compares
 * with the exact certifier's own reading of the same conventions (cell order, local frames).
 *
 *   npx tsx packages/boundary/scripts/gen-containers.ts
 */
import { writeFileSync } from "node:fs";
import { dirname, join } from "node:path";
import { fileURLToPath } from "node:url";
import { POLYTOPES, getPolytope } from "../src/index";

const out: Record<string, unknown> = {};
for (const name of POLYTOPES) {
  const p = getPolytope(name);
  // every cell of the small containers, every 7th of the 600-cell and every 3rd of the 120-cell
  const step = p.cells.length > 200 ? 7 : p.cells.length > 60 ? 3 : 1;
  out[name] = {
    vertices: p.vertices,
    cellCount: p.cells.length,
    cells: p.cells
      .filter((c) => c.index % step === 0)
      .map((c) => ({ index: c.index, ids: c.vertexIds, origin: c.origin, frame: c.frame, normal: c.normal })),
    hingeAngle: p.hingeAngle,
  };
}
const file = join(
  dirname(fileURLToPath(import.meta.url)),
  "..",
  "..",
  "..",
  "python",
  "boundary",
  "fixtures",
  "containers.json",
);
writeFileSync(file, JSON.stringify(out));
console.log("wrote", file);
