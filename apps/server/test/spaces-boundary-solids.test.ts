/**
 * The Platonic solids on the regular 4-polytopes through the server's spaces routes (docs/spaces/BOUNDARY.md 4C): the domain's problems
 * without records and their bounds, submissions of the certified test packings (rotated pieces: q is read, kept and judged), refusals, the strict
 * check a problem without a record needs, and the heavy-check classification by the measured price of a piece.
 */
import { readFileSync } from "node:fs";
import { join, resolve } from "node:path";
import { boundaryIsHeavy, boundaryCheckMs } from "@polypack/spaces";
import { SpaceProblemResponse, SpaceSubmitResponse, SpacesResponse } from "@polypack/spaces";
import { loadSpacesIndex } from "@polypack/spaces/load";
import { beforeEach, describe, expect, it } from "vitest";
import { createApp } from "../src/app";
import { SqliteStore } from "../src/sqlite-store";
import { withUnits } from "./units";

const root = resolve(__dirname, "../../..");
const index = loadSpacesIndex(root);
const packings: { id: string; piece: string; P: number; score_cert: number; config: any }[] = JSON.parse(
  readFileSync(join(root, "python/boundary/fixtures/platonic_packings.json"), "utf8"),
).packings;
const packing = (id: string) => packings.find((p) => p.id === id)!;
const scaled = (c: any, f: number) => ({
  s: c.s * f,
  pieces: c.pieces.map((q: any) => ({ ...q, u: q.u.map((x: number) => x * f) })),
});

let app: ReturnType<typeof createApp>;
let requests = 0;
beforeEach(() => {
  const store = new SqliteStore(":memory:");
  app = createApp({
    store,
    records: [],
    adminToken: "secret",
    clientIp: () => `10.2.${Math.floor(++requests / 200)}.${requests % 200}`,
    log: () => {},
    spaces: { index, store },
  });
});
const get = async (path: string) => {
  const r = await app.request(path);
  return { status: r.status, body: (await r.json()) as any };
};
async function submit(id: string, config: unknown): Promise<{ status: number; body: SpaceSubmitResponse }> {
  const r = await app.request("/api/spaces/submit", {
    method: "POST",
    headers: { "content-type": "application/json" },
    body: JSON.stringify(withUnits({ space: "boundary", id, config, name: "ann" })),
  });
  return { status: r.status, body: (await r.json()) as SpaceSubmitResponse };
}

describe("the solids in the boundary domain", () => {
  it("lists 816 problems on the 4-polytopes (the d-cubes 4,288, the d-simplexes, d-orthoplexes and d-demicubes 432), 32 with records, and their pieces by container", async () => {
    const r = (await get("/api/spaces")).body as SpacesResponse;
    const b = r.spaces.find((s) => s.space === "boundary")!;
    expect(b).toMatchObject({ problems: 816 + 4288 + 432, records: 32 });
    const piece = b.domain!.fields.find((f) => f.key === "piece")!;
    const all = ["ball", "cube", "tetrahedron", "octahedron", "dodecahedron", "icosahedron"];
    for (const n of ["tesseract", "5cell", "16cell", "24cell", "120cell", "600cell"]) expect(piece.dependent!.values![n]).toEqual(all);
    expect(piece.dependent!.values!["7cube"]).toEqual(["ball", "cube"]);
    // the simplex, orthoplex and demicube take balls only
    for (const n of ["5simplex", "10orthoplex", "8demicube"]) expect(piece.dependent!.values![n]).toEqual(["ball"]);
    expect(b.conventions.pose).toMatch(/rotates the standard solid/);
  });

  it("answers a solid's problem with its bound, no record and no configuration, and a rotated cube is a problem of the 24-cell", async () => {
    for (const id of ["boundary/24cell/octahedron/5", "boundary/24cell/cube/3", "boundary/tesseract/icosahedron/16", "boundary/120cell/dodecahedron/1"]) {
      const { status, body } = await get(`/api/spaces/boundary/problems/${id}`);
      expect(status, id).toBe(200);
      const r = body as SpaceProblemResponse;
      expect(r).toMatchObject({ id, status: "no_record", certified: null, certified_config: null, accepts_submissions: true });
      expect(r.bound!.value).toBeGreaterThan(0);
      expect(r.bound!.label).toMatch(/^\[CALC\] (volume|one piece)/);
    }
    // the tesseract's cube keeps its record and its bound
    const cube = (await get("/api/spaces/boundary/problems/boundary/tesseract/cube/9")).body as SpaceProblemResponse;
    expect(cube.status).not.toBe("no_record");
    expect(cube.certified!.score).toBeGreaterThan(0);
    for (const id of ["boundary/6cube/tetrahedron/3", "boundary/24cell/sphere/3"]) expect((await get(`/api/spaces/boundary/problems/${id}`)).status, id).toBe(400);
    expect((await get("/api/spaces/boundary/problems/boundary/24cell/octahedron/25")).status).toBe(404);
  });
});

describe("POST /api/spaces/submit: solids", () => {
  it("accepts every certified test packing, with its rotations, as the first entry, and keeps only a better one after", async () => {
    for (const r of packings) {
      const out = await submit(r.id, r.config);
      expect(out.status, r.id).toBe(200);
      expect(out.body, `${r.id}: ${out.body.reasons?.join("; ")}`).toMatchObject({ valid: true, leader: true, score: r.score_cert });
      const again = await submit(r.id, r.config);
      expect(again.body, r.id).toMatchObject({ valid: true, leader: false, id: null });
    }
    const r = packing("boundary/16cell/cube/2");
    const p = (await get(`/api/spaces/boundary/problems/${r.id}`)).body as SpaceProblemResponse;
    expect(p.site_best!.score).toBe(r.score_cert);
    expect(p.site_config).toEqual(r.config);
  });

  it("refuses a shrunk packing (an edge through a solid, or an overlap), a bound violation, a bad rotation, a rotated ball and a wrong number of pieces", async () => {
    const r = packing("boundary/24cell/tetrahedron/3");
    const shrunk = await submit(r.id, scaled(r.config, 1 - 1e-6));
    expect(shrunk.status).toBe(422);
    expect(shrunk.body.reasons[0]).toMatch(/Tetrahedr(a|on)|edge|overlap/);
    const b = (await get(`/api/spaces/boundary/problems/${r.id}`)).body.bound.value;
    const below = await submit(r.id, { s: b * 0.9, pieces: r.config.pieces });
    expect(below.status).toBe(422);
    expect(below.body.reasons[0]).toBe(`the size ${b * 0.9} is below the lower bound ${b} of this problem`);
    const piece = r.config.pieces[0];
    const rest = r.config.pieces.slice(1);
    expect((await submit(r.id, { s: r.config.s, pieces: [{ ...piece, q: [2, 0, 0, 0] }, ...rest] })).status).toBe(400);
    expect((await submit(r.id, { s: r.config.s, pieces: [{ ...piece, q: [1, 0, 0] }, ...rest] })).status).toBe(400);
    expect((await submit(r.id, { s: r.config.s, pieces: [{ ...piece, cell: 24 }, ...rest] })).status).toBe(400);
    expect((await submit(r.id, { s: r.config.s, pieces: rest })).status).toBe(422);
    // a rotation on a ball is no pose; a ball's id with a solid's configuration
    expect((await submit("boundary/24cell/ball/3", r.config)).status).toBe(400);
    expect((await submit("boundary/24cell/tetrahedron/30", r.config)).status).toBe(400);
  });

  it("needs tolerance 1e-11 of an entry on a problem without a record", async () => {
    // the 24 octahedra of the 24-cell, 5e-10 relative too close: valid at 1e-9, not at 1e-11
    const r = packing("boundary/24cell/octahedron/24");
    const loose = scaled(r.config, (1 - 5e-10) / r.config.s);
    const v = await submit(r.id, loose);
    expect(v.status).toBe(422);
    expect(v.body.reasons[0]).toMatch(/no research record yet/);
  });
});

describe("the heavy-check classification of solids", () => {
  it("prices a piece by its container and queues a check of 50 ms or more", () => {
    const heavy = (polytope: string, piece: string, P: number) => boundaryIsHeavy({ polytope: polytope as never, piece: piece as never, P });
    expect(heavy("120cell", "dodecahedron", 4)).toBe(false);
    expect(heavy("120cell", "dodecahedron", 5)).toBe(true);
    expect(heavy("tesseract", "icosahedron", 12)).toBe(false);
    expect(heavy("tesseract", "icosahedron", 13)).toBe(true);
    expect(heavy("16cell", "cube", 9)).toBe(false);
    expect(heavy("16cell", "cube", 10)).toBe(true);
    expect(heavy("600cell", "tetrahedron", 24)).toBe(false);
    // the aligned cubes and the balls keep their prices
    expect(heavy("tesseract", "cube", 16)).toBe(false);
    expect(boundaryCheckMs({ polytope: "120cell", piece: "dodecahedron", P: 24 })).toBe(264);
  });
});
