import { SPACE_NAMES } from "@polypack/spaces";
import { ND_DIMENSIONS, ndProblemCount } from "@polypack/ndspace/catalogue";
import { boundaryDomain } from "@polypack/boundary";
import { describe, expect, it } from "vitest";
import {
  DIMENSIONS,
  availableFamilies,
  familyFilters,
  familyKey,
  familyOfFilters,
  familyOfSpace,
  familyTitle,
  isAvailable,
  menuGroups,
  scopeHasId,
  whereIs,
} from "../src/dimension-table";
import { BOARD_SPACES } from "../src/space-address";
import { SPACE_FIELDS } from "../src/spaces-model";

/** The labels in the table's order: the four of the first releases, then 4D in 4D, 4D on 5D, 5D in 5D, ... 9D on 10D, 10D in 10D. */
const LABELS = [
  "2D in 2D",
  "2D on 3D",
  "3D in 3D",
  "3D on 4D",
  "4D in 4D",
  "4D on 5D",
  "5D in 5D",
  "5D on 6D",
  "6D in 6D",
  "6D on 7D",
  "7D in 7D",
  "7D on 8D",
  "8D in 8D",
  "8D on 9D",
  "9D in 9D",
  "9D on 10D",
  "10D in 10D",
];

describe("the table of dimensions", () => {
  it("names each dimension by the pieces' dimension, in or on, and the container's", () => {
    expect(DIMENSIONS.map((d) => d.label)).toEqual(LABELS);
    for (const d of DIMENSIONS) {
      const m = /^(\d+)D (in|on) (\d+)D$/.exec(d.label)!;
      expect(m, d.label).not.toBeNull();
      // "in" a container of the same dimension; "on" a boundary of one dimension more
      expect(Number(m[3]) - Number(m[1]), d.label).toBe(m[2] === "in" ? 0 : 1);
      expect(d.id).toBe(`${m[1]}-${m[2]}-${m[3]}`);
      expect(d.families.length, d.label).toBeGreaterThan(0);
    }
  });

  it("files every kind the board can show, Flat first, one family per dimension of the nD spaces", () => {
    const shown = DIMENSIONS.flatMap((d) => availableFamilies(d).map((f) => f.space));
    expect(shown[0]).toBe("polygon");
    // the nD spaces back several families each: ndim seven, boundary three-on-four and the six d-cubes
    expect(new Set(shown)).toEqual(new Set(["polygon", ...BOARD_SPACES]));
    expect(shown.filter((s) => s === "ndim")).toHaveLength(7);
    expect(shown.filter((s) => s === "boundary")).toHaveLength(7);
    expect(whereIs("polygon")).toMatchObject({ dimension: { label: "2D in 2D" }, family: { name: "Flat" } });
    expect(whereIs("curved")?.dimension.label).toBe("2D in 2D");
    for (const s of ["surface", "sphere", "torus"]) expect(whereIs(s)?.dimension.label, s).toBe("2D on 3D");
    expect(whereIs("boundary")).toMatchObject({ dimension: { id: "3-on-4" }, family: { name: "Polytope boundaries" } });
    expect(whereIs("polyhedra")).toMatchObject({ dimension: { label: "3D in 3D" }, family: { name: "Polyhedra" } });
    expect(whereIs("balls")).toMatchObject({ dimension: { label: "3D in 3D" }, family: { name: "Balls" } });
    expect(familyTitle("balls")).toBe("3D in 3D, Balls");
    expect(familyTitle("sphere")).toBe("2D on 3D, Sphere");
    expect(familyTitle("boundary")).toBe("3D on 4D, Polytope boundaries");
  });

  it("files every space of the API, so that the dimensions page can group its records", () => {
    for (const s of SPACE_NAMES) {
      expect(familyOfSpace(s), s).not.toBeNull();
      expect(Object.hasOwn(SPACE_FIELDS, s), s).toBe(true);
    }
  });

  it("lists no dimension as still to come: every one has a kind the board can show", () => {
    const soon = DIMENSIONS.filter((d) => !isAvailable(d));
    expect(soon.map((d) => d.label)).toEqual([]);
    expect(familyOfSpace("polyhedra")?.dimension.label).toBe("3D in 3D");
    expect(familyOfSpace("boundary")?.dimension.label).toBe("3D on 4D");
  });

  it("uses distinct menu keys and distinct names within a dimension, and a short name and a sentence for each family", () => {
    const families = DIMENSIONS.flatMap((d) => d.families);
    expect(new Set(families.map(familyKey)).size).toBe(families.length);
    expect(new Set(DIMENSIONS.map((d) => d.id)).size).toBe(DIMENSIONS.length);
    for (const d of DIMENSIONS) expect(new Set(d.families.map((f) => f.name)).size, d.label).toBe(d.families.length);
    // the menu's text is "<dimension>: <name>", and no two options read alike
    const menu = menuGroups().flatMap((g) => g.options.map((o) => o.text));
    expect(new Set(menu).size).toBe(menu.length);
    for (const f of families) {
      expect(f.short.length, familyKey(f)).toBeGreaterThan(0);
      expect(f.summary, familyKey(f)).toMatch(/\.$/);
    }
    for (const d of DIMENSIONS) expect(d.summary, d.label).toMatch(/\.$/);
  });

  it("gives the board one menu: a group per dimension, an option per family, the ones to come disabled", () => {
    expect(menuGroups().map((g) => g.label)).toEqual(LABELS);
    expect(
      menuGroups()
        .slice(0, 4)
        .flatMap((g) => g.options.map((o) => o.text)),
    ).toEqual([
      "2D in 2D: Flat",
      "2D in 2D: Curved",
      "2D on 3D: Polyhedral surfaces",
      "2D on 3D: Sphere",
      "2D on 3D: Torus",
      "3D in 3D: Polyhedra",
      "3D in 3D: Balls",
      "3D on 4D: Polytope boundaries",
    ]);
    // the nD dimensions, one option each
    expect(menuGroups()[4].options).toEqual([
      { text: "4D in 4D: Polytopes and balls", value: "ndim-4", disabled: false },
    ]);
    expect(menuGroups()[5].options).toEqual([
      { text: "4D on 5D: Polytope boundaries", value: "boundary-5cube", disabled: false },
    ]);
    expect(menuGroups()[16].options.map((o) => o.value)).toEqual(["ndim-10"]);
    // the value is what the board shows: "polygon", the space, or (where a space backs several families) the family's key
    expect(menuGroups()[1].options.map((o) => o.value)).toEqual(["surface", "sphere", "torus"]);
    expect(menuGroups()[2].options.map((o) => o.value)).toEqual(["polyhedra", "balls"]);
    expect(menuGroups()[3].options.map((o) => o.value)).toEqual(["boundary"]);
    expect(
      menuGroups()
        .flatMap((g) => g.options)
        .filter((o) => o.disabled),
    ).toEqual([]);
    const soon = menuGroups([
      {
        id: "4-in-4",
        label: "4D in 4D",
        summary: "Polytopes in a polytope.",
        families: [
          { space: "polytopes", name: "Polytopes", short: "Polytopes", summary: "Soon.", available: false },
          { space: "polygon", name: "Flat", short: "Flat", summary: "Here.", available: true },
        ],
      },
    ]);
    expect(soon).toEqual([
      {
        label: "4D in 4D",
        options: [
          { text: "4D in 4D: Polytopes (coming soon)", value: "polytopes", disabled: true },
          { text: "4D in 4D: Flat", value: "polygon", disabled: false },
        ],
      },
    ]);
  });
});

describe("one space behind several dimensions", () => {
  it("opens the dimension that holds an ndim problem, by its d", () => {
    for (const d of ND_DIMENSIONS) {
      const id = `ndim/${d}/cube/ball/5`;
      expect(whereIs("ndim", id)?.dimension.label, id).toBe(`${d}D in ${d}D`);
      expect(familyKey(whereIs("ndim", id)!.family)).toBe(`ndim-${d}`);
      // the menu's key names it too
      expect(whereIs(`ndim-${d}`)?.dimension.id).toBe(`${d}-in-${d}`);
    }
    // a bare #space=ndim opens 4D in 4D
    expect(whereIs("ndim")?.dimension.label).toBe("4D in 4D");
  });

  it("opens the dimension that holds a boundary problem, by its container", () => {
    // every existing 3D on 4D link keeps working: the tesseract and the five 4-polytopes
    for (const polytope of ["tesseract", "5cell", "16cell", "24cell", "120cell", "600cell"])
      expect(whereIs("boundary", `boundary/${polytope}/ball/3`)?.dimension.label, polytope).toBe("3D on 4D");
    expect(whereIs("boundary")?.dimension.label).toBe("3D on 4D");
    expect(whereIs("boundary", "boundary/tesseract/ball/9")).toMatchObject({ family: { name: "Polytope boundaries" } });
    // the d-cubes, simplexes, orthoplexes and demicubes: "4D on 5D" to "9D on 10D"
    for (let d = 5; d <= 10; d++) {
      const ids = [`boundary/${d}cube/cube/7`, `boundary/${d}simplex/ball/7`, `boundary/${d}orthoplex/ball/7`, `boundary/${d}demicube/ball/7`];
      for (const id of ids) expect(whereIs("boundary", id)?.dimension.label, id).toBe(`${d - 1}D on ${d}D`);
      expect(whereIs(`boundary-${d}cube`)?.dimension.id).toBe(`${d - 1}-on-${d}`);
      expect(whereIs(`boundary-${d}cube`)?.family).toMatchObject({ name: "Polytope boundaries" });
    }
  });

  it("holds every id of both spaces in exactly one family", () => {
    const ids = boundaryDomain.sample(100000);
    expect(ids.length).toBeGreaterThan(4000);
    for (const id of ids) {
      const holders = DIMENSIONS.flatMap((d) => d.families).filter(
        (f) => f.space === "boundary" && scopeHasId(f.scope, id),
      );
      expect(holders, id).toHaveLength(1);
    }
    for (const d of ND_DIMENSIONS) {
      const holders = DIMENSIONS.flatMap((x) => x.families).filter(
        (f) => f.space === "ndim" && scopeHasId(f.scope, `ndim/${d}/ortho/ball/2`),
      );
      expect(holders.map(familyKey), String(d)).toEqual([`ndim-${d}`]);
    }
    // an id outside a family's scope is not held by it
    expect(scopeHasId(DIMENSIONS.find((x) => x.id === "5-in-5")!.families[0].scope, "ndim/6/cube/cube/3")).toBe(false);
  });

  it("opens the records of a dimension from the filters of the table", () => {
    expect(familyFilters(whereIs("ndim-7")!.family)).toEqual({ d: "7" });
    expect(familyFilters(whereIs("boundary-6cube")!.family)).toEqual({ polytope: "6cube" });
    // 3D on 4D holds six polytopes: no single filter opens it
    expect(familyFilters(whereIs("boundary")!.family)).toEqual({});
    expect(familyOfFilters("ndim", { d: "9" })?.dimension.label).toBe("9D in 9D");
    expect(familyOfFilters("ndim", {})?.dimension.label).toBe("4D in 4D");
    expect(familyOfFilters("ndim", { d: "all" })?.dimension.label).toBe("4D in 4D");
    expect(familyOfFilters("boundary", { polytope: "8cube" })?.dimension.label).toBe("7D on 8D");
    expect(familyOfFilters("boundary", { polytope: "24cell" })?.dimension.label).toBe("3D on 4D");
    expect(familyOfFilters("boundary", { polytope: "all" })?.dimension.label).toBe("3D on 4D");
    expect(familyOfFilters("sphere", {})?.dimension.label).toBe("2D on 3D");
    expect(familyOfFilters("planet", {})).toBeNull();
  });

  it("starts each family on a problem of its own scope", () => {
    for (const f of DIMENSIONS.flatMap((d) => d.families)) {
      if (!f.scope) continue;
      expect(f.start, familyKey(f)).toBeDefined();
      expect(scopeHasId(f.scope, f.start!), familyKey(f)).toBe(true);
    }
  });

  it("counts the problems of the nD dimensions as the catalogue does", () => {
    expect(ND_DIMENSIONS.map(ndProblemCount).reduce((a, b) => a + b, 0)).toBe(6131);
  });
});
