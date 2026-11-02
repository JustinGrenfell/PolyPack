/**
 * Balls on the d-simplex, d-orthoplex and d-demicube (d = 5 to 10: 4D on 5D to 9D on 10D) through the server's spaces routes (the real
 * catalogues, an in-memory store): the domain in the listing, problems without records and their bounds, submissions of the certified test
 * packings (python/boundary/fixtures/nd_packings.json), the strict check a problem without a record needs, malformed input, the heavy-check
 * queue and the cap on the work (docs/spaces/BOUNDARY.md 4D).
 */
import { readFileSync } from "node:fs";
import { join, resolve } from "node:path";
import { ND_MAX_WORK, ND_NAMES, ndBound, ndParse } from "@polypack/boundary";
import { boundaryCheckMs, boundaryIsHeavy } from "@polypack/spaces";
import { SpaceProblemResponse, SpaceSubmitResponse, SpacesResponse } from "@polypack/spaces";
import { loadSpacesIndex } from "@polypack/spaces/load";
import { beforeEach, describe, expect, it } from "vitest";
import { createApp } from "../src/app";
import { SqliteStore } from "../src/sqlite-store";
import { withUnits } from "./units";

const root = resolve(__dirname, "../../..");
const index = loadSpacesIndex(root);
interface Packing {
  id: string;
  polytope: string;
  P: number;
  config: { s: number; pieces: { cell: number; u: number[] }[] };
}
const packings: Packing[] = JSON.parse(
  readFileSync(join(root, "python/boundary/fixtures/nd_packings.json"), "utf8"),
).packings;
const packing = (id: string) => packings.find((p) => p.id === id)!;

let app: ReturnType<typeof createApp>;
let requests = 0;
beforeEach(() => {
  const store = new SqliteStore(":memory:");
  app = createApp({
    store,
    records: [],
    adminToken: "secret",
    clientIp: () => `10.3.${Math.floor(++requests / 200)}.${requests % 200}`,
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

describe("the simplex, orthoplex and demicube domain on the server", () => {
  it("is listed among the boundary problems: 18 containers, balls only, P = 1 to 24", async () => {
    const r = (await get("/api/spaces")).body as SpacesResponse;
    const b = r.spaces.find((s) => s.space === "boundary")!;
    const fields = b.domain!.fields;
    for (const n of ND_NAMES) {
      expect(fields[0].values).toContain(n);
      expect(fields[1].dependent!.values![n]).toEqual(["ball"]);
      expect(fields[2].dependent!.ranges![n]).toEqual([1, 24]);
    }
    expect(b.conventions.polytope).toMatch(/5simplex.*10demicube/);
  });

  it("answers a problem without a record with its exact bound and no configuration", async () => {
    for (const [id, name, P] of [
      ["boundary/7simplex/ball/3", "7simplex", 3],
      ["boundary/10orthoplex/ball/24", "10orthoplex", 24],
      ["boundary/8demicube/ball/1", "8demicube", 1],
    ] as const) {
      const { status, body } = await get(`/api/spaces/boundary/problems/${id}`);
      expect(status).toBe(200);
      const r = body as SpaceProblemResponse;
      expect(r).toMatchObject({
        id,
        status: "no_record",
        certified: null,
        certified_config: null,
        accepts_submissions: true,
      });
      const p = ndParse(name)!;
      expect(r.bound!.value).toBe(ndBound(p.family, p.d, P));
      expect(r.bound!.label).toMatch(/^\[CALC\]/);
    }
    const one = (await get("/api/spaces/boundary/problems/boundary/9simplex/ball/1")).body as SpaceProblemResponse;
    // the one-ball bound sqrt(d (d-2)) / 2, as the largest binary64 not above it
    expect(Math.abs(one.bound!.value - Math.sqrt(63) / 2)).toBeLessThan(1e-15);
    expect(one.bound!.value).toBeLessThanOrEqual(Math.sqrt(63) / 2 + 1e-15);
  });

  it("has no problem for an unwritten P, a P above 24, cubes or an unknown container", async () => {
    for (const id of [
      "boundary/7simplex/ball/03",
      "boundary/7simplex/cube/3",
      "boundary/11simplex/ball/3",
      "boundary/7simplex/ball/0",
    ])
      expect((await get(`/api/spaces/boundary/problems/${id}`)).status, id).toBe(400);
    for (const id of ["boundary/7simplex/ball/25", "boundary/10demicube/ball/9999"])
      expect((await get(`/api/spaces/boundary/problems/${id}`)).status, id).toBe(404);
  });
});

describe("POST /api/spaces/submit: the simplex, orthoplex and demicube", () => {
  it("accepts every certified test packing (72), strictly valid at the 1e-11 a problem without a record needs", async () => {
    for (const p of packings) {
      const out = await submit(p.id, p.config);
      expect(out.status, p.id).toBe(200);
      expect(out.body.valid, `${p.id}: ${out.body.reasons.join("; ")}`).toBe(true);
    }
  });

  it("keeps the first valid entry of a problem without a record, then only one that is better", async () => {
    const p = packing("boundary/7demicube/ball/24");
    const first = await submit(p.id, p.config);
    expect(first.body).toMatchObject({ valid: true, leader: true, score: p.config.s });
    const again = await submit(p.id, p.config);
    expect(again.body).toMatchObject({ valid: true, leader: false, id: null });
    const q = (await get(`/api/spaces/boundary/problems/${p.id}`)).body as SpaceProblemResponse;
    expect(q.site_best!.score).toBe(p.config.s);
    expect(q.site_config).toEqual(p.config);
  });

  it("refuses shrunk packings, a bound violation, a bad cell, a rotation, a wrong number of pieces and cubes", async () => {
    const p = packing("boundary/8orthoplex/ball/24");
    const f = 1 - 1e-6;
    const shrunk = await submit(p.id, {
      s: p.config.s * f,
      pieces: p.config.pieces.map((q) => ({ cell: q.cell, u: q.u.map((x) => x * f) })),
    });
    expect(shrunk.status).toBe(422);
    expect(shrunk.body.reasons[0]).toMatch(/overlap|volume|\(d-3\)-face/);
    const b = (await get(`/api/spaces/boundary/problems/${p.id}`)).body.bound.value;
    const below = await submit(p.id, { s: b * 0.9, pieces: p.config.pieces });
    expect(below.status).toBe(422);
    expect(below.body.reasons[0]).toBe(`the size ${b * 0.9} is below the lower bound ${b} of this problem`);
    const rest = p.config.pieces.slice(1);
    expect(
      (await submit(p.id, { s: p.config.s, pieces: [{ ...p.config.pieces[0], cell: 256 }, ...rest] })).status,
    ).toBe(400);
    expect(
      (await submit(p.id, { s: p.config.s, pieces: [{ ...p.config.pieces[0], u: [0, 0, 0] }, ...rest] })).status,
    ).toBe(400);
    const rotated = await submit(p.id, {
      s: p.config.s,
      pieces: [{ ...p.config.pieces[0], q: [1, 0, 0, 0] }, ...rest],
    });
    expect(rotated.status).toBe(400);
    expect(rotated.body.reasons.join(" ")).toMatch(/no rotation q/);
    expect(rotated.body.reasons.join(" ")).toMatch(/is a ball and has no rotation q/);
    expect(rotated.body.reasons.join(" ")).not.toMatch(/axis-aligned/);
    expect((await submit(p.id, { s: p.config.s, pieces: rest })).status).toBe(422);
    expect((await submit("boundary/8orthoplex/ball/25", p.config)).status).toBe(400);
    expect((await submit("boundary/8orthoplex/cube/24", p.config)).status).toBe(400);
  });

  it("needs tolerance 1e-11 of an entry on a problem without a record", async () => {
    // two balls at adjacent cells' centres of the 6-simplex, 5e-10 relative too close: valid at 1e-9, not at 1e-11
    const s0 = Math.sqrt(15);
    const loose = {
      s: s0 * (1 - 5e-10),
      pieces: [
        { cell: 0, u: new Array(5).fill(0) },
        { cell: 1, u: new Array(5).fill(0) },
      ],
    };
    const v = await submit("boundary/6simplex/ball/2", loose);
    expect(v.status).toBe(422);
    expect(v.body.reasons[0]).toMatch(/no research record yet/);
  });
});

describe("the heavy-check classification and the cap on the work", () => {
  it("queues the checks of 50 ms or more by their estimated cost: the 6-demicube from P = 3, the 7-demicube from P = 2, the 8-demicube always, the 10-demicube from P = 10", () => {
    const heavy = (polytope: string, P: number) => boundaryIsHeavy({ polytope: polytope as never, piece: "ball", P });
    expect(heavy("6demicube", 2)).toBe(false);
    expect(heavy("6demicube", 3)).toBe(true);
    expect(heavy("7demicube", 1)).toBe(false);
    expect(heavy("7demicube", 2)).toBe(true);
    expect(heavy("8demicube", 1)).toBe(true);
    expect(heavy("10demicube", 9)).toBe(false);
    expect(heavy("10demicube", 10)).toBe(true);
    for (const n of ND_NAMES.filter((x) => !/^(6|7|8|10)demicube$/.test(x))) expect(heavy(n, 24), n).toBe(false);
    expect(boundaryCheckMs({ polytope: "8demicube", piece: "ball", P: 24 })).toBe(125 * 24);
    expect(ND_MAX_WORK).toBe(60_000);
  });

  it("refuses a crowd that overlaps everywhere at once, within the work cap and the 20 reasons", async () => {
    const pieces = Array.from({ length: 24 }, (_, i) => ({ cell: i % 3, u: new Array(9).fill(0.002 * (i % 5)) }));
    const t0 = performance.now();
    const out = await submit("boundary/10simplex/ball/24", { s: 9, pieces });
    expect(performance.now() - t0).toBeLessThan(5000);
    expect(out.status).toBe(422);
    expect(out.body.reasons.length).toBeLessThanOrEqual(20);
  });
});
