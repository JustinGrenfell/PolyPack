/**
 * Balls and axis-aligned cubes on the d-cubes (d = 5 to 10: 4D on 5D to 9D on 10D) through the server's spaces routes (the real
 * catalogues, an in-memory store): the domain in the listing, problems without records, submissions of the certified test packings
 * (python/boundary/fixtures/ncube_packings.json), the strict check a problem without a record needs, malformed input, the heavy-check
 * queue and the caps on the work.
 */
import { readFileSync } from "node:fs";
import { join, resolve } from "node:path";
import { NCUBE_MAX_P, NCUBE_NAMES, ND_NAMES, ncubeCellLayout } from "@polypack/boundary";
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
  piece: "ball" | "cube";
  P: number;
  kind: string;
  config: { s: number; pieces: { cell: number; u: number[] }[] };
}
const packings: Packing[] = JSON.parse(
  readFileSync(join(root, "python/boundary/fixtures/ncube_packings.json"), "utf8"),
).packings;
const packing = (id: string) => packings.find((p) => p.id === id)!;
const rng = (seed: number) => () => {
  seed = (seed * 1664525 + 1013904223) >>> 0;
  return seed / 2 ** 32;
};

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
    clientIp: () => `10.2.${Math.floor(++requests / 200)}.${requests % 200}`,
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
const layout = (d: number, piece: "ball" | "cube", P: number, seed = 1) => {
  const n = d - 1;
  const L = ncubeCellLayout(d, piece, P, rng(seed));
  return {
    s: L.s * (1 + 1e-9),
    pieces: Array.from({ length: P }, (_, i) => ({
      cell: L.st.cells[i],
      u: Array.from(L.st.U.slice(i * n, (i + 1) * n), (x) => x * (1 + 1e-9)),
    })),
  };
};

describe("the d-cube domain on the server", () => {
  it("is listed among the boundary problems with its containers, pieces and P ranges", async () => {
    const r = (await get("/api/spaces")).body as SpacesResponse;
    const b = r.spaces.find((s) => s.space === "boundary")!;
    expect(b.problems).toBe(816 + 2 * (160 + 384 + 4 * 400) + 18 * 24);
    expect(b.records).toBe(32);
    const fields = b.domain!.fields;
    expect(fields[0].values).toEqual(["tesseract", "5cell", "16cell", "24cell", "120cell", "600cell", ...NCUBE_NAMES, ...ND_NAMES]);
    expect(fields[1].dependent!.values!["7cube"]).toEqual(["ball", "cube"]);
    expect(fields[2].dependent!.ranges!["5cube"]).toEqual([1, 160]);
    expect(fields[2].dependent!.ranges!["6cube"]).toEqual([1, 384]);
    expect(fields[2].dependent!.ranges!["10cube"]).toEqual([1, 400]);
  });

  it("answers a problem without a record with its bound and no configuration", async () => {
    const { status, body } = await get("/api/spaces/boundary/problems/boundary/5cube/cube/160");
    expect(status).toBe(200);
    const r = body as SpaceProblemResponse;
    expect(r).toMatchObject({
      id: "boundary/5cube/cube/160",
      status: "no_record",
      certified: null,
      certified_config: null,
      accepts_submissions: true,
    });
    expect(r.fields).toEqual({ polytope: "5cube", piece: "cube", P: 160 });
    expect(r.bound!.value).toBe(2);
    expect(r.bound!.label).toMatch(/^\[CALC\].*volume bound/);
    const one = (await get("/api/spaces/boundary/problems/boundary/10cube/ball/1")).body as SpaceProblemResponse;
    expect(one.bound!.value).toBe(0.7071067811865475);
  });

  it("has no problem for an unwritten P, a P above the container's largest, an unknown container or piece", async () => {
    for (const id of [
      "boundary/5cube/ball/05",
      "boundary/11cube/ball/3",
      "boundary/5cube/sphere/3",
      "boundary/5cube/ball/0",
    ])
      expect((await get(`/api/spaces/boundary/problems/${id}`)).status, id).toBe(400);
    for (const id of [
      "boundary/5cube/ball/161",
      "boundary/6cube/cube/385",
      "boundary/10cube/ball/401",
      "boundary/10cube/ball/9999",
    ])
      expect((await get(`/api/spaces/boundary/problems/${id}`)).status, id).toBe(404);
  });
});

describe("POST /api/spaces/submit: the d-cubes", () => {
  it("accepts every certified test packing (42), strictly valid at the 1e-11 a problem without a record needs", async () => {
    for (const p of packings) {
      const out = await submit(p.id, p.config);
      expect(out.status, p.id).toBe(200);
      expect(out.body.valid, `${p.id}: ${out.body.reasons.join("; ")}`).toBe(true);
    }
  });

  it("keeps the first valid entry of a problem without a record, then only one that is better", async () => {
    const p = packing("boundary/8cube/cube/16");
    const first = await submit(p.id, p.config);
    expect(first.body).toMatchObject({ valid: true, leader: true, score: p.config.s });
    const again = await submit(p.id, p.config);
    expect(again.body).toMatchObject({ valid: true, leader: false, id: null });
    const q = (await get(`/api/spaces/boundary/problems/${p.id}`)).body as SpaceProblemResponse;
    expect(q.site_best!.score).toBe(p.config.s);
    expect(q.site_config).toEqual(p.config);
  });

  it("refuses shrunk packings, a bound violation, a bad cell, a rotation and a wrong number of pieces", async () => {
    const p = packing("boundary/6cube/ball/12");
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
    expect((await submit(p.id, { s: p.config.s, pieces: [{ ...p.config.pieces[0], cell: 12 }, ...rest] })).status).toBe(
      400,
    );
    expect(
      (await submit(p.id, { s: p.config.s, pieces: [{ ...p.config.pieces[0], u: [0, 0, 0] }, ...rest] })).status,
    ).toBe(400);
    const rotated = await submit(p.id, {
      s: p.config.s,
      pieces: [{ ...p.config.pieces[0], q: [1, 0, 0, 0] }, ...rest],
    });
    expect(rotated.status).toBe(400);
    expect(rotated.body.reasons.join(" ")).toMatch(/no rotation q/);
    expect((await submit(p.id, { s: p.config.s, pieces: rest })).status).toBe(422);
    expect((await submit("boundary/6cube/ball/385", p.config)).status).toBe(400);
  });

  it("needs tolerance 1e-11 of an entry on a problem without a record", async () => {
    // two balls at the centres of adjacent cells of the 5-cube, 5e-10 relative too close: valid at 1e-9, not at 1e-11
    const loose = {
      s: 1 - 5e-10,
      pieces: [
        { cell: 0, u: [0, 0, 0, 0] },
        { cell: 2, u: [0, 0, 0, 0] },
      ],
    };
    const v = await submit("boundary/5cube/ball/2", loose);
    expect(v.status).toBe(422);
    expect(v.body.reasons[0]).toMatch(/no research record yet/);
  });

  it("accepts a tiling by cubes at the largest P of the 5-cube and the 6-cube", async () => {
    for (const [d, P] of [
      [5, 160],
      [6, 384],
    ] as const) {
      const n = d - 1;
      const pts = Array.from({ length: 2 ** n }, (_, idx) =>
        Array.from({ length: n }, (_, a) => (((idx >> a) & 1) * 2 - 1) / 2),
      );
      const lam = 1 + 1e-9;
      const pieces = Array.from({ length: 2 * d }, (_, c) =>
        pts.map((u) => ({ cell: c, u: u.map((x) => x * lam) })),
      ).flat();
      expect(pieces).toHaveLength(P);
      const out = await submit(`boundary/${d}cube/cube/${P}`, { s: 2 * lam, pieces });
      expect(out.status, `${d}cube`).toBe(200);
      expect(out.body.valid, out.body.reasons.join("; ")).toBe(true);
    }
  });
});

describe("the heavy-check classification and the caps", () => {
  it("queues balls by their estimated cost: 0.3 to 1.5 ms per ball, 0.02 ms per cube", () => {
    const heavy = (polytope: string, piece: "ball" | "cube", P: number) =>
      boundaryIsHeavy({ polytope: polytope as never, piece, P });
    expect(heavy("5cube", "ball", 160)).toBe(false);
    expect(heavy("5cube", "ball", 167)).toBe(true);
    expect(heavy("10cube", "ball", 33)).toBe(false);
    expect(heavy("10cube", "ball", 34)).toBe(true);
    for (const n of NCUBE_NAMES) expect(heavy(n, "cube", 400), n).toBe(false);
    expect(boundaryCheckMs({ polytope: "10cube", piece: "ball", P: 400 })).toBeLessThan(700);
    expect(NCUBE_MAX_P["10cube"]).toBe(400);
  });

  it("checks 400 balls on the 10-cube through the queue and charges the work", async () => {
    const config = layout(10, "ball", 400);
    const out = await submit("boundary/10cube/ball/400", config);
    expect(out.status).toBe(200);
    expect(out.body.valid, out.body.reasons.join("; ")).toBe(true);
    expect(checks.length).toBeGreaterThan(0);
  });

  it("refuses a crowd that overlaps everywhere at once, within the work cap and the 20 reasons", async () => {
    // 400 balls at the same few points: every pair is a violation; the check stops after 20 reasons
    const pieces = Array.from({ length: 400 }, (_, i) => ({ cell: i % 3, u: new Array(9).fill(0.01 * (i % 5)) }));
    const t0 = performance.now();
    const out = await submit("boundary/10cube/ball/400", { s: 2, pieces });
    expect(performance.now() - t0).toBeLessThan(5000);
    expect(out.status).toBe(422);
    expect(out.body.reasons.length).toBeLessThanOrEqual(20);
  });
});
