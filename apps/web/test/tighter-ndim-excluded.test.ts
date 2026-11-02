/**
 * The nD in nD space (4D to 10D, nothing to do with the game's catalogue): Tighter takes its problems from explicit slices
 * and id patterns, not from the spaces' domains, so the new space never reaches the game. This pins that down.
 */
import { readFileSync, readdirSync } from "node:fs";
import { join } from "node:path";
import { describe, expect, it } from "vitest";
import { SPACE_DOMAINS } from "@polypack/spaces";
import { DIMENSIONS } from "../src/dimension-table";
import { sliceOfId } from "../src/tighter/catalogue";
import { spaceFields } from "../src/tighter/families";
import { familyDimension } from "../src/tighter/frontier";
import { SLICES } from "../src/tighter/slices";

const TIGHTER = join(__dirname, "../src/tighter");

describe("Tighter and the nD spaces", () => {
  it("has no slice for them and does not read their ids as problems", () => {
    expect(Object.keys(SPACE_DOMAINS)).toContain("ndim");
    expect(SLICES.some((s) => /ndim/.test(JSON.stringify(s)))).toBe(false);
    for (const id of ["ndim/4/cube/cube/16", "ndim/8/421/ball/5", "ndim/10/ball/ball/20"]) {
      expect(spaceFields(id), id).toBeNull();
      expect(sliceOfId(id), id).toBeNull();
    }
  });

  it("reads none of the 6,131 ids of the nD in nD space as a problem of the game", () => {
    const ids = SPACE_DOMAINS.ndim!.sample(100000);
    expect(ids).toHaveLength(6131);
    for (const id of ids) {
      expect(spaceFields(id), id).toBeNull();
      expect(sliceOfId(id), id).toBeNull();
    }
  });

  it("keeps the game's Dimensions to the four of the first releases: none of 4D in 4D to 10D in 10D, none of the d-cubes", () => {
    const game = new Set(
      (["polygon", "curved", "surface", "sphere", "torus", "polyhedra", "balls"] as const).map(familyDimension),
    );
    expect([...game].sort()).toEqual(["2-in-2", "2-on-3", "3-in-3"]);
    const site = DIMENSIONS.map((d) => d.id);
    const added = site.filter((id) => !["2-in-2", "2-on-3", "3-in-3", "3-on-4"].includes(id));
    expect(added).toHaveLength(13);
    for (const id of added) expect(game.has(id as never), id).toBe(false);
  });

  it("is not changed by the site's table of dimensions: no file of the game names the nD spaces or the table", () => {
    // the game takes its problems from its own slices and id patterns; a mention of the table, the board's spaces or the
    // nD ids in its source would be a way in
    for (const file of readdirSync(TIGHTER, { recursive: true, encoding: "utf8" }).filter((f) => /\.(ts|css)$/.test(f))) {
      const text = readFileSync(join(TIGHTER, file), "utf8")
        // comments may name them
        .replace(/\/\*[\s\S]*?\*\//g, "")
        .replace(/^\s*\/\/.*$/gm, "");
      expect(text, file).not.toMatch(/dimension-table|space-address|BOARD_SPACES|["'`]ndim[\/"'`]|\bndim\//);
      // the Dimensions it names are the four of the first releases
      for (const m of text.matchAll(/["'`](\d+-(?:in|on)-\d+)["'`]/g))
        expect(["2-in-2", "2-on-3", "3-in-3", "3-on-4"], `${file}: ${m[1]}`).toContain(m[1]);
    }
    // and the page of the game is as it was: it loads no code of the nD views
    const page = readFileSync(join(__dirname, "../tighter.html"), "utf8");
    expect(page).not.toMatch(/ndim|ncube|boundary|space\/render/);
  });
});
