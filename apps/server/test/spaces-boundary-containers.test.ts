/**
 * Balls on the 5-cell, 16-cell, 24-cell, 120-cell and 600-cell through the server's spaces routes (the real catalogues, an
 * in-memory store): the domain in the listing, problems without records, submissions of the certified test packings, the
 * strict check a problem without a record needs, malformed input, and the heavy-check classification.
 */
import { readFileSync } from "node:fs";
import { join, resolve } from "node:path";
import { getPolytope } from "@polypack/boundary";
import { boundaryCheckMs, boundaryIsHeavy } from "@polypack/spaces";
import { SpaceProblemResponse, SpaceSubmitResponse, SpacesResponse } from "@polypack/spaces";
import { loadSpacesIndex } from "@polypack/spaces/load";
import { beforeEach, describe, expect, it } from "vitest";
import { createApp } from "../src/app";
import { SqliteStore } from "../src/sqlite-store";
import { withUnits } from "./units";

const root = resolve(__dirname, "../../..");
const index = loadSpacesIndex(root);
const packings: { id: string; polytope: string; P: number; score_cert: number; config: any }[] = JSON.parse(
  readFileSync(join(root, "python/boundary/fixtures/general_packings.json"), "utf8"),
).packings;
const packing = (id: string) => packings.find((p) => p.id === id)!;
const scaled = (c: any, f: number) => ({
  s: c.s * f,
  pieces: c.pieces.map((q: any) => ({ ...q, u: q.u.map((x: number) => x * f) })),
});

let app: ReturnType<typeof createApp>;
let requests = 0;
let checks: { workMs: number }[] = [];
beforeEach(() => {
  checks = [];
  const store = new SqliteStore(":memory:");
  app = createApp({
    store,
    records: [],
    adminToken: "secret",
    clientIp: () => `10.1.${Math.floor(++requests / 200)}.${requests % 200}`,
    log: () => {},
    spaces: { index, store, onCheck: (c) => checks.push(c) },
  });
});
const get = async (path: string) => {
  const r = await app.request(path);
  return { status: r.status, body: (await r.json()) as any };
};
async function submit(
  id: string,
  config: unknown,
  extra: object = {},
): Promise<{ status: number; body: SpaceSubmitResponse }> {
  const r = await app.request("/api/spaces/submit", {
    method: "POST",
    headers: { "content-type": "application/json" },
    body: JSON.stringify(withUnits({ space: "boundary", id, config, name: "ann", ...extra })),
  });
  return { status: r.status, body: (await r.json()) as SpaceSubmitResponse };
}

describe("the boundary domain on the server", () => {
  it("is listed with 5,104 problems (816 on the 4-polytopes, 4,288 on the d-cubes), 32 of them with records, and its dependent fields", async () => {
    const r = (await get("/api/spaces")).body as SpacesResponse;
    const b = r.spaces.find((s) => s.space === "boundary")!;
    expect(b).toMatchObject({ problems: 816 + 2 * (160 + 384 + 4 * 400), records: 32 });
    expect(b.domain!.fields.map((f) => f.key)).toEqual(["polytope", "piece", "P"]);
    expect(b.domain!.fields[0].values).toEqual([
      "tesseract",
      "5cell",
      "16cell",
      "24cell",
      "120cell",
      "600cell",
      "5cube",
      "6cube",
      "7cube",
      "8cube",
      "9cube",
      "10cube",
    ]);
    expect(b.conventions.polytope).toMatch(/211\.6.*352\.6/);
  });

  it("answers a problem without a record with its bound and no configuration", async () => {
    const { status, body } = await get("/api/spaces/boundary/problems/boundary/24cell/ball/5");
    expect(status).toBe(200);
    const r = body as SpaceProblemResponse;
    expect(r).toMatchObject({
      id: "boundary/24cell/ball/5",
      status: "no_record",
      certified: null,
      certified_config: null,
      accepts_submissions: true,
    });
    expect(r.fields).toEqual({ polytope: "24cell", piece: "ball", P: 5 });
    expect(r.bound!.value).toBeGreaterThan(0.99);
    expect(r.bound!.label).toMatch(/^\[CALC\]/);
  });

  it("has no problem for a solid on a d-cube, an unwritten P or a P above the container's largest", async () => {
    for (const id of ["boundary/6cube/tetrahedron/3", "boundary/24cell/ball/03", "boundary/9cell/ball/3"])
      expect((await get(`/api/spaces/boundary/problems/${id}`)).status, id).toBe(400);
    for (const id of ["boundary/24cell/ball/25", "boundary/tesseract/ball/17", "boundary/600cell/ball/100"])
      expect((await get(`/api/spaces/boundary/problems/${id}`)).status, id).toBe(404);
  });
});

describe("POST /api/spaces/submit: balls on the other containers", () => {
  it("keeps the first valid entry of a problem without a record, then only one that is better", async () => {
    const r = packing("boundary/24cell/ball/6");
    const first = await submit(r.id, r.config);
    expect(first.status).toBe(200);
    expect(first.body).toMatchObject({ valid: true, leader: true, score: r.score_cert });
    expect(first.body.id).not.toBeNull();
    const again = await submit(r.id, r.config);
    expect(again.body).toMatchObject({ valid: true, leader: false, id: null });
    const worse = await submit(r.id, scaled(r.config, 1.01));
    expect(worse.body).toMatchObject({ valid: true, leader: false, id: null });
    const p = (await get(`/api/spaces/boundary/problems/${r.id}`)).body as SpaceProblemResponse;
    expect(p.site_best!.score).toBe(r.score_cert);
    expect(p.site_config).toEqual(r.config);
  });

  it("accepts a certified test packing of every container", async () => {
    for (const r of packings.filter((q) => q.P === 2 || q.P === 3)) {
      const out = await submit(r.id, r.config);
      expect(out.status, r.id).toBe(200);
      expect(out.body.valid, `${r.id}: ${out.body.reasons.join("; ")}`).toBe(true);
    }
  });

  it("refuses shrunk packings, a bound violation, a bad cell, a rotated ball and a wrong number of pieces", async () => {
    const r = packing("boundary/16cell/ball/4");
    const shrunk = await submit(r.id, scaled(r.config, 1 - 1e-6));
    expect(shrunk.status).toBe(422);
    expect(shrunk.body.reasons[0]).toMatch(/overlap|edge/);
    const b = (await get(`/api/spaces/boundary/problems/${r.id}`)).body.bound.value;
    const below = await submit(r.id, { s: b * 0.9, pieces: r.config.pieces });
    expect(below.status).toBe(422);
    expect(below.body.reasons[0]).toBe(`the size ${b * 0.9} is below the lower bound ${b} of this problem`);
    const piece = r.config.pieces[0];
    expect(
      (await submit(r.id, { s: r.config.s, pieces: [{ ...piece, cell: 16 }, ...r.config.pieces.slice(1)] })).status,
    ).toBe(400);
    expect(
      (await submit(r.id, { s: r.config.s, pieces: [{ ...piece, q: [1, 0, 0, 0] }, ...r.config.pieces.slice(1)] }))
        .status,
    ).toBe(400);
    expect((await submit(r.id, { s: r.config.s, pieces: r.config.pieces.slice(1) })).status).toBe(422);
    expect((await submit("boundary/16cell/sphere/4", r.config)).status).toBe(400);
    expect((await submit("boundary/6cube/icosahedron/4", r.config)).status).toBe(400);
    expect((await submit("boundary/16cell/ball/30", r.config)).status).toBe(400);
  });

  it("needs tolerance 1e-11 of an entry on a problem without a record", async () => {
    // two balls at the centres of adjacent 24-cell cells, 5e-10 relative too close: valid at 1e-9, not at 1e-11
    const s0 = Math.sqrt(6) / 2;
    const near = { index: getPolytope("24cell").cells[0].facets[0].neighbor };
    const loose = {
      s: s0 * (1 - 5e-10),
      pieces: [
        { cell: 0, u: [0, 0, 0] },
        { cell: near.index, u: [0, 0, 0] },
      ],
    };
    const v = await submit("boundary/24cell/ball/2", loose);
    expect(v.status).toBe(422);
    expect(v.body.valid).toBe(false);
    expect(v.body.reasons[0]).toMatch(/no research record yet/);
  });
});

describe("the heavy-check classification", () => {
  it("queues balls on the 120-cell from P = 5 and on the 24-cell from P = 22, and nothing else", () => {
    const heavy = (polytope: string, P: number) => boundaryIsHeavy({ polytope: polytope as never, piece: "ball", P });
    expect(heavy("120cell", 4)).toBe(false);
    expect(heavy("120cell", 5)).toBe(true);
    expect(heavy("24cell", 21)).toBe(false);
    expect(heavy("24cell", 22)).toBe(true);
    for (const n of ["tesseract", "5cell", "16cell", "600cell"]) expect(heavy(n, 24), n).toBe(false);
    expect(boundaryCheckMs({ polytope: "120cell", piece: "ball", P: 24 })).toBeLessThan(400);
  });

  it("checks a P = 6 packing on the 120-cell through the queue and charges its work", async () => {
    // six balls at the centres of six cells of the 120-cell, from s0 = 0.4491 (adjacent centres are 2 r s apart)
    const s = 0.4491;
    const out = await submit("boundary/120cell/ball/6", {
      s,
      pieces: Array.from({ length: 6 }, (_, i) => ({ cell: 20 * i, u: [0, 0, 0] })),
    });
    expect(out.status).toBe(200);
    expect(out.body.valid, out.body.reasons.join("; ")).toBe(true);
    expect(checks.length).toBeGreaterThan(0);
  });
});
