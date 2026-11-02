/** The board's picker for the boundary domain: the pieces and the P range depend on the polytope. */
import { boundaryDomain } from "@polypack/boundary";
import { describe, expect, it } from "vitest";
import { BOUNDARY_ROUND, boundaryRound, defaultBudget, freshSize, roundSize } from "../src/space/model";
import { pickFromDomain } from "../src/space/model";

const pick = (wish: Record<string, string>) => pickFromDomain("boundary", boundaryDomain.info, wish);
const values = (p: ReturnType<typeof pick>, key: string) =>
  p.fields.find((f) => f.key === key)!.options.map((o) => o.value);

describe("pickFromDomain on the boundary domain", () => {
  it("opens the tesseract with balls and the Platonic solids, P = 1 to 16, and the default problem", () => {
    const p = pick({});
    expect(p.id).toBe("boundary/tesseract/ball/9");
    expect(values(p, "polytope")).toEqual([
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
    expect(values(p, "piece")).toEqual(["ball", "cube", "tetrahedron", "octahedron", "dodecahedron", "icosahedron"]);
    expect(values(p, "P")).toHaveLength(16);
    expect(p.fields.find((f) => f.key === "polytope")!.options.map((o) => o.text)).toEqual([
      "Tesseract",
      "5-cell",
      "16-cell",
      "24-cell",
      "120-cell",
      "600-cell",
      "5-cube",
      "6-cube",
      "7-cube",
      "8-cube",
      "9-cube",
      "10-cube",
    ]);
    expect(p.fields.find((f) => f.key === "piece")!.options.map((o) => o.text)).toEqual([
      "Balls",
      "Cubes",
      "Tetrahedra",
      "Octahedra",
      "Dodecahedra",
      "Icosahedra",
    ]);
  });

  it("offers balls and cubes on the d-cubes, P = 1 to 160, 384 and 400, and moves a wish to the nearest offered value", () => {
    const p = pick({ polytope: "5cube", piece: "cube", P: "200" });
    expect(p.id).toBe("boundary/5cube/cube/160");
    expect(values(p, "piece")).toEqual(["ball", "cube"]);
    expect(values(p, "P")).toHaveLength(160);
    expect(pick({ polytope: "6cube", piece: "ball", P: "384" }).id).toBe("boundary/6cube/ball/384");
    expect(values(pick({ polytope: "10cube", piece: "ball" }), "P")).toHaveLength(400);
  });

  it("offers balls and the five solids with P = 1 to 24 on the other polytopes, and moves a wish to the nearest offered value", () => {
    const p = pick({ polytope: "24cell", piece: "cube", P: "20" });
    expect(p.id).toBe("boundary/24cell/cube/20");
    expect(values(p, "piece")).toEqual(["ball", "cube", "tetrahedron", "octahedron", "dodecahedron", "icosahedron"]);
    expect(values(p, "P")).toHaveLength(24);
    // a solid is kept when the polytope changes (all six take them); a d-cube takes balls and cubes only
    expect(pick({ polytope: "120cell", piece: "icosahedron", P: "7" }).id).toBe("boundary/120cell/icosahedron/7");
    expect(pick({ polytope: "6cube", piece: "icosahedron", P: "7" }).id).toMatch(/^boundary\/6cube\/(ball|cube)\/7$/);
    expect(pick({ polytope: "120cell", piece: "ball", P: "24" }).id).toBe("boundary/120cell/ball/24");
    // from the 24-cell's P = 24 back to the tesseract: the nearest P it has
    expect(pick({ polytope: "tesseract", piece: "cube", P: "24" }).id).toBe("boundary/tesseract/cube/16");
  });
});

describe("the rounds of 3D on 4D", () => {
  const problem = (polytope: string, P: number) =>
    ({ space: "boundary", id: `boundary/${polytope}/ball/${P}`, P, polytope, piece: "ball" }) as never;

  it("keep the tesseract's, and size the other containers' to about two seconds of the chain metric", () => {
    expect(boundaryRound("tesseract", 16)).toBe(BOUNDARY_ROUND);
    expect(roundSize(problem("tesseract", 9))).toBe(BOUNDARY_ROUND);
    expect(freshSize(problem("tesseract", 9))).toBe(BOUNDARY_ROUND);
    expect(defaultBudget(problem("tesseract", 9))).toEqual({ time_s: 30, evals: 3000000 });
    // 0.8 ms a ball on the 120-cell: P = 8 is 6.4 ms an evaluation, 300 of them 2 s
    expect(boundaryRound("120cell", 8)).toBe(300);
    expect(boundaryRound("24cell", 8)).toBe(1000);
    expect(roundSize(problem("24cell", 8))).toBe(1000);
    expect(freshSize(problem("120cell", 8))).toBe(300);
    expect(defaultBudget(problem("5cell", 3))).toEqual({ time_s: 60, evals: 10 * boundaryRound("5cell", 3) });
    // never fewer than 200, never more than the tesseract's
    expect(boundaryRound("120cell", 64)).toBe(200);
    expect(boundaryRound("5cell", 1)).toBeLessThanOrEqual(BOUNDARY_ROUND);
  });

  it("size the Platonic solids' rounds from their own cost per piece, on the tesseract too", () => {
    const solid = (polytope: string, piece: string, P: number) =>
      ({ space: "boundary", id: `boundary/${polytope}/${piece}/${P}`, P, polytope, piece }) as never;
    // 0.7 ms a piece on the tesseract: P = 4 is 2.8 ms an evaluation, 700 of them 2 s; the aligned cubes keep the tesseract's round
    expect(boundaryRound("tesseract", 4, "icosahedron")).toBe(700);
    expect(boundaryRound("tesseract", 4, "cube")).toBe(BOUNDARY_ROUND);
    // 4 ms a piece on the 120-cell: P = 8 is 32 ms an evaluation, 60 of them 2 s, so the least, 200
    expect(boundaryRound("120cell", 8, "dodecahedron")).toBe(200);
    expect(roundSize(solid("24cell", "octahedron", 8))).toBe(boundaryRound("24cell", 8, "octahedron"));
    expect(defaultBudget(solid("tesseract", "tetrahedron", 3))).toEqual({ time_s: 60, evals: 10 * boundaryRound("tesseract", 3, "tetrahedron") });
    expect(defaultBudget(solid("tesseract", "cube", 3))).toEqual({ time_s: 30, evals: 3000000 });
  });
});
