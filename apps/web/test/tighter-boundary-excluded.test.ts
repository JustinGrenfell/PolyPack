/**
 * The 3D on 4D space (balls on the tesseract and, since 2026-11-01, on the 5-cell, 16-cell, 24-cell, 120-cell and 600-cell; since
 * 2026-11-02 also the five Platonic solids, freely rotated, on all six) and
 * the (d-1)D on dD spaces of the d-cubes (4D on 5D to 9D on 10D, balls and axis-aligned cubes, since 2026-11-01) are
 * not in Tighter: the game takes its problems from explicit slices and id patterns (src/tighter/families.ts spaceFields,
 * slices.ts SLICES), none of which names `boundary`, and not from the spaces' domains. The tesseract was never dealt, and the new
 * containers and pieces cannot be: this pins that down, with no line of the game changed.
 *
 * The expected id count of the 4-polytopes is 816 since 2026-11-02 (it was 152): 664 more problems, the tetrahedron, octahedron,
 * dodecahedron and icosahedron on the tesseract (4 x 16) and the five solids on each of the other five (5 x 5 x 24).
 */
import { describe, expect, it } from "vitest";
import { NCUBE_NAMES, boundaryDomain } from "@polypack/boundary";
import { SPACE_DOMAINS } from "@polypack/spaces";
import { DIMENSIONS } from "../src/dimension-table";
import { sliceOfId } from "../src/tighter/catalogue";
import { spaceFields } from "../src/tighter/families";
import { familyDimension } from "../src/tighter/frontier";
import { SLICES } from "../src/tighter/slices";

describe("Tighter and 3D on 4D", () => {
  it("has no slice for it and does not read any of its 816 ids as a problem", () => {
    expect(SPACE_DOMAINS.boundary).toBe(boundaryDomain);
    expect(SLICES.some((s) => /boundary|tesseract|5cell|16cell|24cell|120cell|600cell/.test(JSON.stringify(s)))).toBe(
      false,
    );
    const ids = boundaryDomain.sample(100000).filter((id) => !/^boundary\/\d+cube\//.test(id));
    expect(ids).toHaveLength(816);
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

describe("Tighter and the boundary's new Dimensions (3D on 4D, 4D on 5D to 9D on 10D)", () => {
  it("deals no problem of the boundary space whatever the Dimension an id opens on the site", () => {
    // the site files each id of the space under a Dimension (dimension-table.ts); the game has no family for any of them
    const ids = boundaryDomain.sample(100000);
    expect(ids).toHaveLength(816 + 4288);
    for (const id of ids) {
      const where = DIMENSIONS.flatMap((d) => d.families.map((f) => ({ d, f }))).find(
        ({ f }) => f.space === "boundary" && f.scope!.values.includes(id.split("/")[1]),
      );
      expect(where, id).toBeDefined();
      expect(spaceFields(id), id).toBeNull();
      expect(sliceOfId(id), id).toBeNull();
    }
    // the game's families map to the first three Dimensions only: 3D on 4D and the d-cubes are not among them
    const game = new Set(
      (["polygon", "curved", "surface", "sphere", "torus", "polyhedra", "balls"] as const).map(familyDimension),
    );
    for (const id of ["3-on-4", "4-on-5", "5-on-6", "6-on-7", "7-on-8", "8-on-9", "9-on-10"]) expect(game.has(id as never), id).toBe(false);
  });
});
