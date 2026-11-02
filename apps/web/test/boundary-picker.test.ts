/** The board's picker for the boundary domain: the pieces and the P range depend on the polytope. */
import { boundaryDomain } from "@polypack/boundary";
import { describe, expect, it } from "vitest";
import { BOUNDARY_ROUND, boundaryRound, defaultBudget, freshSize, roundSize } from "../src/space/model";
import { pickFromDomain } from "../src/space/model";

const pick = (wish: Record<string, string>) => pickFromDomain("boundary", boundaryDomain.info, wish);
const values = (p: ReturnType<typeof pick>, key: string) =>
  p.fields.find((f) => f.key === key)!.options.map((o) => o.value);

describe("pickFromDomain on the boundary domain", () => {
  it("opens the tesseract with balls and cubes, P = 1 to 16, and the default problem", () => {
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
    expect(values(p, "piece")).toEqual(["ball", "cube"]);
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
    expect(p.fields.find((f) => f.key === "piece")!.options.map((o) => o.text)).toEqual(["Balls", "Cubes"]);
  });

  it("offers balls and cubes on the d-cubes, P = 1 to 160, 384 and 400, and moves a wish to the nearest offered value", () => {
    const p = pick({ polytope: "5cube", piece: "cube", P: "200" });
    expect(p.id).toBe("boundary/5cube/cube/160");
    expect(values(p, "piece")).toEqual(["ball", "cube"]);
    expect(values(p, "P")).toHaveLength(160);
    expect(pick({ polytope: "6cube", piece: "ball", P: "384" }).id).toBe("boundary/6cube/ball/384");
    expect(values(pick({ polytope: "10cube", piece: "ball" }), "P")).toHaveLength(400);
  });

  it("offers balls only and P = 1 to 24 on the other polytopes, and moves a wish to the nearest offered value", () => {
    const p = pick({ polytope: "24cell", piece: "cube", P: "20" });
    expect(p.id).toBe("boundary/24cell/ball/20");
    expect(values(p, "piece")).toEqual(["ball"]);
    expect(values(p, "P")).toHaveLength(24);
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
});
