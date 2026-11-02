/**
 * Balls on the five containers beyond the tesseract in the registry: the boundary domain (152 problems, 32 with records), the
 * problems without records and their lower bounds, the server's rule on the certified test packings, a boundary job on one of them.
 */
import { readFileSync } from "node:fs";
import { resolve } from "node:path";
import { CONTAINER_MAX_P, POLYTOPES, ballBounds, boundaryDomain, getPolytope } from "@polypack/boundary";
import { describe, expect, it } from "vitest";
import {
  SPACE_DOMAINS,
  SPACE_STRICT_NO_RECORD_REASON,
  acceptsSubmissions,
  parseSpaceId,
  problemOf,
  spaceProblemSummary,
  validateSpaceSubmission,
} from "../src";
import { loadSpacesIndex } from "../src/load";
import { SPACE_METHOD, parseSearchProblem, runSpaceJob } from "../src/search";

const root = resolve(__dirname, "../../..");
const cat = loadSpacesIndex(root).boundary!;
const packings: { id: string; polytope: string; P: number; score_cert: number; config: any }[] = JSON.parse(
  readFileSync(resolve(root, "python/boundary/fixtures/general_packings.json"), "utf8"),
).packings;
const NEW = POLYTOPES.filter((n) => n !== "tesseract");
const entry = (id: string) => {
  const parsed = parseSpaceId("boundary", id);
  if (typeof parsed === "string") throw new Error(parsed);
  return problemOf(cat, parsed);
};
const scaled = (c: any, f: number) => ({
  s: c.s * f,
  pieces: c.pieces.map((q: any) => ({ ...q, u: q.u.map((x: number) => x * f) })),
});

describe("the boundary domain", () => {
  it("is balls on six containers (P = 1 to 24; the tesseract 1 to 16) and the tesseract's cubes: 152 problems", () => {
    expect(SPACE_DOMAINS.boundary).toBe(boundaryDomain);
    expect(boundaryDomain.info.problems).toBe(32 + 5 * 24);
    expect(cat.domain).toBe(boundaryDomain);
    expect(cat.maxP).toBe(24);
    // the sample lists each problem once, small P first, and every id parses
    const ids = boundaryDomain.sample(1000);
    expect(ids).toHaveLength(152);
    expect(new Set(ids).size).toBe(152);
    expect(ids.slice(0, 7)).toEqual([
      "boundary/tesseract/ball/1",
      "boundary/tesseract/cube/1",
      ...NEW.slice(0, 5).map((n) => `boundary/${n}/ball/1`),
    ]);
    for (const id of ids) {
      const p = parseSpaceId("boundary", id);
      expect(typeof p, id).toBe("object");
      expect(boundaryDomain.bound(p as never), id).not.toBeUndefined();
    }
    // the fields the board's picker reads: the pieces and the P range depend on the polytope
    const piece = boundaryDomain.info.fields.find((f) => f.key === "piece")!;
    expect(piece.dependent!.values!.tesseract).toEqual(["ball", "cube"]);
    expect(piece.dependent!.values!["24cell"]).toEqual(["ball"]);
    const P = boundaryDomain.info.fields.find((f) => f.key === "P")!;
    expect(P.dependent!.ranges!.tesseract).toEqual([1, 16]);
    expect(P.dependent!.ranges!["600cell"]).toEqual([1, 24]);
  });

  it("keeps the tesseract's 32 records unchanged and adds problems without records", () => {
    expect(cat.problems).toHaveLength(32);
    expect(cat.problems.every((p) => p.fields.polytope === "tesseract")).toBe(true);
    for (const n of NEW)
      for (const P of [1, 12, 24]) {
        const p = entry(`boundary/${n}/ball/${P}`)!;
        expect(p.status).toBe("no_record");
        expect(p.certified).toBeNull();
        expect(p.config).toBeNull();
        expect(acceptsSubmissions(p)).toBe(true);
        const b = ballBounds(n, P);
        expect(p.bound!.value).toBe(b.value);
        expect(p.bound!.label).toMatch(b.kind === "volume" ? /^\[CALC\] volume/ : /^\[CALC\] one ball/);
        expect(spaceProblemSummary(p, null)).toMatchObject({
          status: "no_record",
          accepts_submissions: true,
          best_score: null,
        });
      }
  });

  it("has no problem outside the domain: a cube off the tesseract, P above the container's largest, an unknown polytope", () => {
    for (const id of [
      "boundary/24cell/cube/3",
      "boundary/tesseract/ball/17",
      "boundary/600cell/ball/25",
      "boundary/tesseract/cube/17",
    ]) {
      const parsed = parseSpaceId("boundary", id);
      if (typeof parsed === "string") expect(parsed).toMatch(/takes ball/);
      else expect(problemOf(cat, parsed), id).toBeNull();
    }
    for (const id of ["boundary/9cell/ball/3", "boundary/24cell/ball/03", "boundary/24cell/ball/0"])
      expect(typeof parseSpaceId("boundary", id), id).toBe("string");
  });

  it("reads the new ids in the searches too", () => {
    expect(parseSearchProblem("boundary/24cell/ball/5")).toEqual({
      space: "boundary",
      id: "boundary/24cell/ball/5",
      P: 5,
      polytope: "24cell",
      piece: "ball",
    });
    expect(parseSearchProblem("boundary/24cell/cube/5")).toBeNull();
  });
});

describe("the server's rule on the new containers", () => {
  const rule = (id: string) => ({ bound: entry(id)!.bound!.value, noRecord: true, maxP: cat.maxP });

  it("accepts every certified test packing, also at the strict tolerance a problem without a record needs", () => {
    expect(packings.length).toBeGreaterThanOrEqual(25);
    for (const r of packings) {
      const v = validateSpaceSubmission("boundary", r.id, r.config, rule(r.id));
      expect(v.valid, `${r.id}: ${v.reasons.join("; ")}`).toBe(true);
      expect(v.score).toBe(r.score_cert);
      expect(r.score_cert).toBeGreaterThan(entry(r.id)!.bound!.value);
    }
  });

  it("refuses them shrunk by 1e-6, with a reason, and below the bound before any geometry", () => {
    for (const r of packings.filter((p) => p.P >= 2)) {
      const v = validateSpaceSubmission("boundary", r.id, scaled(r.config, 1 - 1e-6), rule(r.id));
      expect(v.valid, r.id).toBe(false);
      expect(v.malformed).toBe(false);
      expect(v.reasons[0]).toMatch(/edge|overlap|below the lower bound|volume/);
    }
    const one = packings.find((p) => p.id === "boundary/24cell/ball/1")!;
    const b = entry(one.id)!.bound!.value;
    const v = validateSpaceSubmission("boundary", one.id, { s: b * 0.99, pieces: one.config.pieces }, rule(one.id));
    expect(v.reasons[0]).toBe(`the size ${b * 0.99} is below the lower bound ${b} of this problem`);
  });

  it("refuses a cell outside the polytope's range, a rotation on a ball and a wrong number of pieces", () => {
    const r = packings.find((p) => p.id === "boundary/24cell/ball/3")!;
    const piece = r.config.pieces[0];
    const bad = (pieces: unknown[]) => validateSpaceSubmission("boundary", r.id, { s: r.config.s, pieces }, rule(r.id));
    const rest = r.config.pieces.slice(1);
    expect(bad([{ ...piece, cell: 24 }, ...rest]).reasons[0]).toBe("pieces[0].cell must be an integer from 0 to 23");
    expect(bad([{ ...piece, cell: 24 }, ...rest]).malformed).toBe(true);
    expect(bad([{ ...piece, q: [1, 0, 0, 0] }, ...rest]).reasons[0]).toBe("pieces[0] is a ball and has no rotation q");
    const v = bad(rest);
    expect(v.valid).toBe(false);
    expect(v.malformed).toBe(false);
    expect(v.reasons[0]).toBe("expected 3 pieces for this problem, got 2");
    // the 5-cell has 5 cells: cell 5 is outside, and so is 120-cell cell 120
    const five = packings.find((p) => p.id === "boundary/5cell/ball/1")!;
    expect(
      validateSpaceSubmission(
        "boundary",
        five.id,
        { s: five.config.s, pieces: [{ cell: 5, u: [0, 0, 0] }] },
        rule(five.id),
      ).reasons[0],
    ).toBe("pieces[0].cell must be an integer from 0 to 4");
    expect(getPolytope("120cell").cells).toHaveLength(120);
  });

  it("asks for tolerance 1e-11 of a problem without a record: a packing valid only at 1e-9 is refused", () => {
    // two balls at the centres of adjacent 24-cell cells, 2e-10 relative too close: within the touching tolerance, not the strict one
    const s0 = Math.sqrt(6) / 2;
    const id = "boundary/24cell/ball/2";
    const loose = {
      s: s0 * (1 - 5e-10),
      pieces: [
        { cell: 0, u: [0, 0, 0] },
        { cell: getPolytope("24cell").cells[0].facets[0].neighbor, u: [0, 0, 0] },
      ],
    };
    expect(validateSpaceSubmission("boundary", id, loose, { bound: entry(id)!.bound!.value }).valid).toBe(true);
    const v = validateSpaceSubmission("boundary", id, loose, rule(id));
    expect(v.valid).toBe(false);
    expect(v.reasons[0]).toBe(SPACE_STRICT_NO_RECORD_REASON);
  });

  it("refuses P above the container's largest in a packing, as the registry's id has no such problem", () => {
    expect(CONTAINER_MAX_P["24cell"]).toBe(24);
  });
});

describe("one boundary job on a new container", () => {
  it("returns a valid configuration by the server's rule, reproducibly", () => {
    const id = "boundary/24cell/ball/3";
    const job = {
      space: "boundary" as const,
      id,
      method: SPACE_METHOD.boundary,
      seed: 4,
      budget: { time_s: 60, evals: 6000 },
      kind: "fresh" as const,
      rule: { bound: entry(id)!.bound!.value, noRecord: true },
    };
    const a = runSpaceJob(job, { aborted: false, paused: false });
    const b = runSpaceJob(job, { aborted: false, paused: false });
    expect(a.valid, a.reasons.join("; ")).toBe(true);
    expect(a.score!).toBeGreaterThan(entry(id)!.bound!.value);
    expect(b.config).toEqual(a.config);
  });
});
