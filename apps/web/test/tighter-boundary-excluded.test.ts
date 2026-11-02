/**
 * The 3D on 4D space (balls on the tesseract and, since 2026-11-01, on the 5-cell, 16-cell, 24-cell, 120-cell and 600-cell) is
 * not in Tighter: the game takes its problems from explicit slices and id patterns (src/tighter/families.ts spaceFields,
 * slices.ts SLICES), none of which names `boundary`, and not from the spaces' domains. The tesseract was never dealt, and the new
 * containers cannot be: this pins that down, with no line of the game changed.
 */
import { describe, expect, it } from "vitest";
import { boundaryDomain } from "@polypack/boundary";
import { SPACE_DOMAINS } from "@polypack/spaces";
import { sliceOfId } from "../src/tighter/catalogue";
import { spaceFields } from "../src/tighter/families";
import { SLICES } from "../src/tighter/slices";

describe("Tighter and 3D on 4D", () => {
  it("has no slice for it and does not read any of its 152 ids as a problem", () => {
    expect(SPACE_DOMAINS.boundary).toBe(boundaryDomain);
    expect(SLICES.some((s) => /boundary|tesseract|5cell|16cell|24cell|120cell|600cell/.test(JSON.stringify(s)))).toBe(
      false,
    );
    const ids = boundaryDomain.sample(1000);
    expect(ids).toHaveLength(152);
    for (const id of ids) {
      expect(spaceFields(id), id).toBeNull();
      expect(sliceOfId(id), id).toBeNull();
    }
  });
});
