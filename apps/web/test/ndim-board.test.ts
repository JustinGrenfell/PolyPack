/**
 * The pure parts of the nD dimensions on the site: the board's picker of one dimension, the Records page's tabs, pickers and
 * counts of one dimension, the d-cube configuration reader, and the addresses.
 */
import { SPACE_DOMAINS } from "@polypack/spaces";
import { describe, expect, it } from "vitest";
import { DIMENSIONS, familyKey, whereIs } from "../src/dimension-table";
import {
  domainPRange,
  kindHref,
  pickers,
  recordsSource,
  scopedDomainCount,
  seriesTitle,
  tabModel,
} from "../src/records-model";
import { parseSpaceAddress, spaceAddress } from "../src/space-address";
import { parseSearchProblem } from "@polypack/spaces/search";
import { isQueuedCheck, pickFromDomain, submitAnswer, submittingText } from "../src/space/model";
import { facetName, ncubeDimension, readNcubeConfig, worldCentre } from "../src/space/render/ncube";
import { SPACE_FIELDS } from "../src/spaces-model";

const ndim = SPACE_DOMAINS.ndim!.info;
const boundary = SPACE_DOMAINS.boundary!.info;
const scopeOf = (key: string) => whereIs(key)!.family.scope!;
const info = (space: "ndim" | "boundary") =>
  ({ space, domain: space === "ndim" ? ndim : boundary }) as unknown as Parameters<typeof scopedDomainCount>[0];

describe("the board's picker of one dimension", () => {
  it("offers the dimension's shapes and range of P, and fixes d", () => {
    const pick = pickFromDomain("ndim", ndim, { d: "6", piece: "221", container: "122", P: "30" }, scopeOf("ndim-6"));
    expect(pick.id).toBe("ndim/6/221/122/24");
    const d = pick.fields.find((f) => f.key === "d")!;
    expect(d.options.map((o) => o.value)).toEqual(["6"]);
    expect(d.hidden).toBe(true);
    expect(pick.fields.filter((f) => f.hidden).map((f) => f.key)).toEqual(["d"]);
    // the shapes of 6D: its four regular ones, the Gosset polytopes and the ball
    expect(pick.fields.find((f) => f.key === "piece")!.options.map((o) => o.value)).toEqual([
      "simplex",
      "cube",
      "ortho",
      "demi",
      "221",
      "122",
      "ball",
    ]);
    // titles are words
    expect(pick.fields.find((f) => f.key === "piece")!.options[1].text).toBe("Hypercube");
    expect(pick.fields.find((f) => f.key === "P")!.options).toHaveLength(23);
  });

  it("opens a problem of another dimension at the nearest values the dimension has", () => {
    // 24cell is a 4D shape: 5D keeps the starting hypercube; P 36 is above 5D's range of 32
    const pick = pickFromDomain("ndim", ndim, { d: "4", piece: "24cell", container: "600cell", P: "36" }, scopeOf("ndim-5"));
    expect(pick.id).toBe("ndim/5/cube/cube/32");
    const four = pickFromDomain("ndim", ndim, { d: "4", piece: "24cell", container: "600cell", P: "36" }, scopeOf("ndim-4"));
    expect(four.id).toBe("ndim/4/24cell/600cell/36");
    // no scope: as before, every dimension
    expect(pickFromDomain("ndim", ndim, { d: "9", piece: "ball", container: "ball", P: "7" }).fields[0].options).toHaveLength(7);
    for (const d of [4, 5, 6, 7, 8, 9, 10])
      expect(pickFromDomain("ndim", ndim, {}, scopeOf(`ndim-${d}`)).id.startsWith(`ndim/${d}/`), String(d)).toBe(true);
  });

  it("limits 3D on 4D to its six polytopes and a d-cube to itself", () => {
    const three = pickFromDomain("boundary", boundary, { polytope: "5cube" }, scopeOf("boundary"));
    expect(three.fields.find((f) => f.key === "polytope")!.options.map((o) => o.value)).toEqual([
      "tesseract",
      "5cell",
      "16cell",
      "24cell",
      "120cell",
      "600cell",
    ]);
    expect(three.fields.find((f) => f.key === "polytope")!.hidden).toBeUndefined();
    // the tesseract with cubes has 16 pieces at most; a link with other fields keeps the polytope it names
    expect(pickFromDomain("boundary", boundary, { polytope: "24cell", piece: "ball", P: "7" }, scopeOf("boundary")).id).toBe(
      "boundary/24cell/ball/7",
    );
    for (let d = 5; d <= 10; d++) {
      const pick = pickFromDomain("boundary", boundary, { polytope: "tesseract", piece: "cube", P: "9" }, scopeOf(`boundary-${d}cube`));
      expect(pick.id, String(d)).toBe(`boundary/${d}cube/cube/9`);
      expect(pick.fields.find((f) => f.key === "polytope")!.hidden).toBe(true);
    }
  });
});

describe("the Records page of one dimension", () => {
  it("links a dimension to its filter, and marks the current one", () => {
    expect(kindHref("ndim", { d: "5" })).toBe("/records.html#space=ndim&d=5");
    expect(kindHref("torus")).toBe("/records.html#space=torus");
    const { dimensions, kinds } = tabModel("ndim", { d: "6" });
    expect(dimensions.map((d) => d.text)).toEqual(DIMENSIONS.map((d) => d.label));
    expect(dimensions.filter((d) => d.current).map((d) => d.text)).toEqual(["6D in 6D"]);
    expect(dimensions.find((d) => d.text === "4D on 5D")!.href).toBe("/records.html#space=boundary&polytope=5cube");
    expect(dimensions.find((d) => d.text === "3D on 4D")!.href).toBe("/records.html#space=boundary");
    expect(dimensions.find((d) => d.text === "10D in 10D")!.href).toBe("/records.html#space=ndim&d=10");
    // one kind in each: no second row of tabs
    expect(kinds).toEqual([]);
    // the current tab follows the d-cube, and an old link (#space=boundary&polytope=tesseract) stays on 3D on 4D
    expect(tabModel("boundary", { polytope: "7cube" }).dimensions.filter((d) => d.current).map((d) => d.text)).toEqual(["6D on 7D"]);
    expect(tabModel("boundary", { polytope: "tesseract" }).dimensions.filter((d) => d.current).map((d) => d.text)).toEqual(["3D on 4D"]);
    expect(tabModel("boundary").dimensions.filter((d) => d.current).map((d) => d.text)).toEqual(["3D on 4D"]);
    expect(tabModel("ndim").dimensions.filter((d) => d.current).map((d) => d.text)).toEqual(["4D in 4D"]);
  });

  it("offers a dimension's own values in the pickers", () => {
    const four = pickers("ndim", [], info("ndim"), scopeOf("ndim-4"));
    expect(four.find((p) => p.spec.key === "d")!.options.map((o) => o.value)).toEqual(["4"]);
    expect(four.find((p) => p.spec.key === "piece")!.options.map((o) => o.text)).toEqual(
      expect.arrayContaining(["Hypercube", "Orthoplex", "24-cell", "600-cell", "Ball"]),
    );
    const cube = pickers("boundary", [], info("boundary"), scopeOf("boundary-8cube"));
    expect(cube.find((p) => p.spec.key === "polytope")!.options).toEqual([{ value: "8cube", text: "8-cube" }]);
    const four4 = pickers("boundary", [], info("boundary"), scopeOf("boundary"));
    expect(four4.find((p) => p.spec.key === "polytope")!.options.map((o) => o.value)).toHaveLength(6);
  });

  it("counts a dimension's problems from the domain: 6,131 in all, 152 and 4,288 on the boundaries", () => {
    const per = [4, 5, 6, 7, 8, 9, 10].map((d) => scopedDomainCount(info("ndim"), scopeOf(`ndim-${d}`)));
    expect(per).toEqual([1911, 775, 1127, 684, 684, 475, 475]);
    expect(per.reduce((a, b) => a! + b!, 0)).toBe(6131);
    expect(scopedDomainCount(info("ndim"), null)).toBe(6131);
    expect(scopedDomainCount(info("boundary"), scopeOf("boundary"))).toBe(152);
    const cubes = [5, 6, 7, 8, 9, 10].map((d) => scopedDomainCount(info("boundary"), scopeOf(`boundary-${d}cube`)));
    expect(cubes).toEqual([320, 768, 800, 800, 800, 800]);
    expect(cubes.reduce((a, b) => a! + b!, 0)).toBe(4288);
    expect(scopedDomainCount(undefined, null)).toBeNull();
  });

  it("knows the range of P of the chosen dimension", () => {
    expect(domainPRange(info("ndim"), { d: "4" })).toEqual([2, 40]);
    expect(domainPRange(info("ndim"), { d: "5" })).toEqual([2, 32]);
    expect(domainPRange(info("ndim"), { d: "10" })).toEqual([2, 20]);
    expect(domainPRange(info("boundary"), { polytope: "5cube" })).toEqual([1, 160]);
  });

  it("titles a family as a sentence and says where the records come from", () => {
    const specs = SPACE_FIELDS.ndim;
    expect(seriesTitle(specs, { d: "4", piece: "cube", container: "ortho" })).toBe("4D: Hypercubes in an orthoplex");
    expect(seriesTitle(specs, { d: "6", piece: "ball", container: "221" })).toBe("6D: Balls in a 2_21 polytope");
    expect(seriesTitle(specs, { d: "8", piece: "24cell", container: "ball" })).toBe("8D: 24-cells in a ball");
    expect(recordsSource("ndim", { d: "5" })).toMatch(/Zarzuelo.*certified them again/);
    expect(recordsSource("ndim", { d: "5" })).toMatch(/Nakajima/);
    expect(recordsSource("ndim", { d: "9" })).toMatch(/^No research record yet/);
    expect(recordsSource("sphere", {})).toBeNull();
  });
});

describe("addresses", () => {
  it("reads and writes an ndim problem, and keeps the boundary ones", () => {
    expect(parseSpaceAddress("#space=ndim/5/cube/ball/3")).toEqual({ space: "ndim", id: "ndim/5/cube/ball/3" });
    expect(parseSpaceAddress("#space=ndim")).toEqual({ space: "ndim", id: null });
    expect(spaceAddress("ndim", "ndim/9/ortho/ortho/4")).toBe("#space=ndim/9/ortho/ortho/4");
    for (const id of ["boundary/tesseract/ball/9", "boundary/600cell/ball/2", "boundary/6cube/cube/12"])
      expect(parseSpaceAddress(`#space=${id}`)).toEqual({ space: "boundary", id });
    // the dimension an address opens
    expect(familyKey(whereIs("ndim", "ndim/9/ortho/ortho/4")!.family)).toBe("ndim-9");
    expect(familyKey(whereIs("boundary", "boundary/6cube/cube/12")!.family)).toBe("boundary-6cube");
    expect(familyKey(whereIs("boundary", "boundary/600cell/ball/2")!.family)).toBe("boundary");
  });
});

describe("the d-cube's configuration", () => {
  it("names facets, places centres in R^d and refuses malformed configurations", () => {
    expect(ncubeDimension("7cube")).toBe(7);
    expect(ncubeDimension("tesseract")).toBeNull();
    // cell 2i + (sigma > 0): the facet x_i = sigma s/2
    expect(facetName(0)).toBe("−x1");
    expect(facetName(1)).toBe("+x1");
    expect(facetName(5)).toBe("+x3");
    // the facet's axis is at +-s/2 and the other coordinates are u in order
    expect(worldCentre({ cell: 3, u: [0.1, 0.2, 0.3, 0.4] }, 5, 2)).toEqual([0.1, 1, 0.2, 0.3, 0.4]);
    expect(worldCentre({ cell: 0, u: [0.1, 0.2, 0.3, 0.4] }, 5, 2)).toEqual([-1, 0.1, 0.2, 0.3, 0.4]);
    const ok = { s: 2, pieces: [{ cell: 9, u: [0, 0, 0, 0] }] };
    expect(readNcubeConfig(ok, 5)).toEqual(ok);
    for (const bad of [
      null,
      { s: 2 },
      { s: -1, pieces: [] },
      { s: 2, pieces: [{ cell: 10, u: [0, 0, 0, 0] }] },
      { s: 2, pieces: [{ cell: 0, u: [0, 0, 0] }] },
      { s: 2, pieces: [{ cell: 0, u: [0, 0, 0, NaN] }] },
      { s: 2, pieces: [[0, 0, 0, 0]] },
    ])
      expect(readNcubeConfig(bad, 5), JSON.stringify(bad)).toBeNull();
    // a tesseract configuration is not a 5-cube's
    expect(readNcubeConfig({ s: 2, pieces: [{ cell: 0, u: [0, 0, 0] }] }, 5)).toBeNull();
  });
});

describe("an entry whose check the server queues", () => {
  it("says it may wait, and shows the server's answer", () => {
    const small = parseSearchProblem("ndim/4/cube/cube/3")!;
    const big = parseSearchProblem("ndim/8/421/421/20")!;
    expect(isQueuedCheck(small)).toBe(false);
    expect(isQueuedCheck(big)).toBe(true);
    expect(submittingText(small)).toBe("Submitting…");
    expect(submittingText(big)).toMatch(/one at a time.*wait for its turn/);
    // the answers of a queued check: kept, valid but not kept, a full queue, a timeout
    expect(submitAnswer(200, { valid: true, score: 3.5, leader: true, id: 7 }, "ndim")).toMatchObject({ kind: "kept" });
    expect(submitAnswer(200, { valid: true, score: 3.5, leader: false, id: null, best_score: 3.4 }, "ndim").text).toMatch(
      /^✓ Valid, score 3\.500000, but not kept/,
    );
    const busy = submitAnswer(503, { reasons: ["the server is busy checking other entries; try again in a few seconds"] }, "ndim");
    expect(busy.kind).toBe("busy");
    expect(busy.text).toBe("✗ The server is busy checking other entries; try again in a few seconds. Nothing was stored; try again later.");
    expect(submitAnswer(422, { reasons: ["the check took more than 8 s of processor time"] }, "ndim")).toMatchObject({ kind: "invalid" });
  });
});
