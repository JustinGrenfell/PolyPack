/**
 * The Platonic solids through the spaces registry (docs/spaces/BOUNDARY.md 4C): one boundary job on a solid, the generic methods on the solid kernel, the
 * registry's own check of a solid's configuration against the kernel's, and the listing of the methods a solid problem offers.
 */
import { readFileSync } from "node:fs";
import { resolve } from "node:path";
import { listMethods } from "@polypack/methods";
import { describe, expect, it } from "vitest";
import { type SpaceJob, SPACE_METHOD, checkSpaceConfig, parseSearchProblem, runSpaceJob, spaceKernel, spaceMethodsFor } from "../src/search";
import { validateSpaceSubmission } from "../src/validate";
import { AUTO } from "../src/search";

const packings: { id: string; score_cert: number; config: any }[] = JSON.parse(
  readFileSync(resolve(__dirname, "../../../python/boundary/fixtures/platonic_packings.json"), "utf8"),
).packings;
const packing = (id: string) => packings.find((p) => p.id === id)!;
const go = { aborted: false, paused: false };
const job = (id: string, o: Partial<SpaceJob> = {}): SpaceJob => ({
  space: "boundary",
  id,
  method: SPACE_METHOD.boundary,
  seed: 3,
  budget: { time_s: 60, evals: 3000 },
  kind: "fresh",
  rule: {},
  ...o,
});

describe("one boundary job on a solid", () => {
  it("finds a valid packing from nothing, by the server's rule, reproducibly", () => {
    for (const id of ["boundary/tesseract/octahedron/1", "boundary/5cell/cube/1"]) {
      const a = runSpaceJob(job(id), go);
      expect(a.valid, `${id}: ${a.reasons.join("; ")}`).toBe(true);
      expect(checkSpaceConfig("boundary", id, a.config, {})).toMatchObject({ valid: true, score: a.score });
      const b = runSpaceJob(job(id), go);
      expect(b.config).toEqual(a.config);
    }
  }, 120000);

  it("improves an enlarged certified test packing", () => {
    const r0 = packing("boundary/16cell/cube/2");
    const f = 1.2;
    const start = { s: r0.config.s * f, pieces: r0.config.pieces.map((p: any) => ({ ...p, u: p.u.map((x: number) => x * f) })) };
    const out = runSpaceJob(job("boundary/16cell/cube/2", { kind: "iterate", start, budget: { time_s: 60, evals: 4000 } }), go);
    expect(out.valid, out.reasons.join("; ")).toBe(true);
    expect(out.score!).toBeLessThan(start.s);
  }, 120000);
});

describe("the registry's rule and the kernel's agree on a solid", () => {
  it("accepts the certified packings, refuses a shrunk one, and gives the same score", () => {
    for (const id of ["boundary/tesseract/tetrahedron/3", "boundary/24cell/tetrahedron/3", "boundary/120cell/dodecahedron/2"]) {
      const r = packing(id);
      const reg = validateSpaceSubmission("boundary", id, r.config, { bound: null, certified: null });
      expect(reg.valid, `${id}: ${reg.reasons.join("; ")}`).toBe(true);
      const k = spaceKernel(parseSearchProblem(id)!)!;
      const kv = k.validate(r.config);
      expect(kv.valid, kv.reasons.join("; ")).toBe(true);
      expect(kv.score).toBe(reg.score);
      const shrunk = { s: r.config.s * (1 - 1e-6), pieces: r.config.pieces.map((p: any) => ({ ...p, u: p.u.map((x: number) => x * (1 - 1e-6)) })) };
      const regBad = validateSpaceSubmission("boundary", id, shrunk, { bound: null, certified: null });
      const kBad = k.validate(shrunk);
      expect(regBad.valid, id).toBe(false);
      expect(kBad.valid, id).toBe(false);
      expect(kBad.reasons).toEqual(regBad.reasons);
    }
  });
});

describe("the generic methods on the solids' kernel", () => {
  it("lists the methods its capabilities allow (seeds, swap; no soft shapes)", () => {
    for (const id of ["boundary/24cell/octahedron/4", "boundary/tesseract/icosahedron/3", "boundary/5cell/cube/5"]) {
      const p = parseSearchProblem(id)!;
      const k = spaceKernel(p)!;
      expect(k.dim).toBe(6);
      expect(k.capabilities).toEqual(["seeds", "swap"]);
      const expected = listMethods().filter((m) => m.requires.every((c) => k.capabilities.includes(c))).map((m) => m.name);
      expect(spaceMethodsFor(p).map((m) => m.name), id).toEqual([AUTO, ...expected]);
      expect(expected).toEqual(expect.arrayContaining(["penalty", "mc_compress", "asc", "basin_hopping", "shrink_repair", "seed"]));
      expect(expected).not.toContain("harden");
    }
  });

  it("run end to end: penalty, mc_compress and basin_hopping each return a valid packing or none", () => {
    const id = "boundary/tesseract/tetrahedron/1";
    for (const method of ["penalty", "mc_compress", "basin_hopping"]) {
      const r = runSpaceJob(job(id, { method, budget: { time_s: 20, evals: 1500 }, params: {} }), go);
      if (r.config === null) continue;
      expect(r.valid, `${method}: ${r.reasons.join("; ")}`).toBe(true);
      expect(checkSpaceConfig("boundary", id, r.config, {})).toMatchObject({ valid: true });
    }
  }, 120000);
});
