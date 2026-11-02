/**
 * The 3D on 4D space (balls on the tesseract and, since 2026-11-01, on the 5-cell, 16-cell, 24-cell, 120-cell and 600-cell) and
 * the (d-1)D on dD spaces of the d-cubes (4D on 5D to 9D on 10D, balls and axis-aligned cubes, since 2026-11-01) are
 * not in Tighter: the game takes its problems from explicit slices and id patterns (src/tighter/families.ts spaceFields,
 * slices.ts SLICES), none of which names `boundary`, and not from the spaces' domains. The tesseract was never dealt, and the new
 * containers cannot be: this pins that down, with no line of the game changed.
 */
import { describe, expect, it } from "vitest";
import { NCUBE_NAMES, boundaryDomain } from "@polypack/boundary";
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
    const ids = boundaryDomain.sample(100000).filter((id) => !/^boundary\/\d+cube\//.test(id));
    expect(ids).toHaveLength(152);
    for (const id of ids) {
      expect(spaceFields(id), id).toBeNull();
      expect(sliceOfId(id), id).toBeNull();
    }
  });
});

describe("Tighter and (d-1)D on dD, d = 5 to 10", () => {
  it("has no slice for the d-cubes and does not read any of their 4,288 ids as a problem", () => {
    expect(
      SLICES.some(
        (s) => /(^|[^0-9])(5|6|7|8|9|10)cube/.test(JSON.stringify(s)) || /boundary\/\d+cube/.test(JSON.stringify(s)),
      ),
    ).toBe(false);
    const ids = boundaryDomain.sample(100000).filter((id) => /^boundary\/\d+cube\//.test(id));
    // balls and cubes for P up to 160, 384 and 400 (four times)
    expect(ids).toHaveLength(2 * (160 + 384 + 4 * 400));
    for (const n of NCUBE_NAMES) expect(ids.some((id) => id.startsWith(`boundary/${n}/ball/`))).toBe(true);
    for (const id of ids) {
      expect(spaceFields(id), id).toBeNull();
      expect(sliceOfId(id), id).toBeNull();
    }
  });
});
