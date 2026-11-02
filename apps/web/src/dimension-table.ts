/**
 * The dimensions of PolyPack's problems and the families filed under each: one table, read by the
 * board's Problem set menu (dimension-switch.ts) and by the dimensions page.
 *
 * A dimension is named by the pieces' dimension, "in" or "on", and the container's dimension.
 * "in" means inside a container; "on" means on a closed surface or boundary. Polygons in polygons
 * are 2D in 2D; polygons on the surface of a solid are 2D on 3D.
 *
 * Each family names the board mode that shows it: "polygon" (the classic board) or a space of the
 * API (/api/spaces/<space>). The API calls a dimension's families "spaces", and the space names
 * stay as they are, so that links (#space=...) and downloaded volunteer clients keep working.
 *
 * To make a family available, set `available: true`; its space must then be one the board can
 * show (BOARD_SPACES in space-address.ts), which the type below checks. A family that is not
 * available is listed, disabled, as "coming soon".
 *
 * No imports at run time, so that the default board can read it without loading the spaces.
 */
import type { BoardSpace } from "./space-address";

/** What a family opens on the board: the classic polygon board, or a space mode. */
export type FamilyTarget = BoardSpace | "polygon";

/**
 * One API space can back several dimensions (ndim: 4D in 4D to 10D in 10D; boundary: 3D on 4D and the d-cubes). A family
 * with a scope shows only the problems whose id has one of `values` at position `at` ("ndim/5/..." has "5" at 1;
 * "boundary/5cube/..." has "5cube" at 1). A family without one shows all of its space.
 */
export interface FamilyScope {
  /** The field's name in the records table's filters ("d", "polytope"). */
  field: string;
  /** Its position in the problem id. */
  at: number;
  values: readonly string[];
}

interface FamilyText {
  /**
   * What the Problem set menu's option carries, and what names the family in code; the space when absent. Needed where
   * one space backs several families ("ndim-5").
   */
  key?: string;
  /** Which problems of the space this family shows; all of them when absent. */
  scope?: FamilyScope;
  /** The problem the board opens the family with when a link names none (a family with a scope; else the board's default of the space). */
  start?: string;
  /** The name in the board's Problem set menu and on the dimensions page. */
  name: string;
  /** A shorter name, for the records table's picker. */
  short: string;
  /** One plain sentence. */
  summary: string;
}

export type Family =
  | (FamilyText & { space: FamilyTarget; available: true })
  /** Coming soon: `space` is the name the API will use. */
  | (FamilyText & { space: string; available: false });

export interface Dimension {
  /** The id, e.g. "2-on-3": the board's <main> carries it as data-dimension. */
  id: string;
  /** "2D on 3D". */
  label: string;
  /** One or two plain sentences. */
  summary: string;
  families: readonly Family[];
}


// ---------- the nD dimensions (4D in 4D to 10D in 10D, 4D on 5D to 9D on 10D) ----------

/** The shapes of each dimension of the ndim space, in words (docs/spaces/NDIM.md section 2), and its largest P. */
const ND_TEXT: Record<number, { shapes: string; maxP: number }> = {
  4: { shapes: "the 5-cell, tesseract, 16-cell, 24-cell, 120-cell and 600-cell", maxP: 40 },
  5: { shapes: "the 5-simplex, 5-cube, 5-orthoplex and 5-demicube", maxP: 32 },
  6: { shapes: "the 6-simplex, 6-cube, 6-orthoplex, 6-demicube and the Gosset polytopes 2_21 and 1_22", maxP: 24 },
  7: { shapes: "the 7-simplex, 7-cube, 7-orthoplex, 7-demicube and the Gosset polytope 3_21", maxP: 20 },
  8: { shapes: "the 8-simplex, 8-cube, 8-orthoplex, 8-demicube and the Gosset polytope 4_21", maxP: 20 },
  9: { shapes: "the 9-simplex, 9-cube, 9-orthoplex and 9-demicube", maxP: 20 },
  10: { shapes: "the 10-simplex, 10-cube, 10-orthoplex and 10-demicube", maxP: 20 },
};

/** "4D in 4D" and so on: equal polytopes or balls in a polytope or ball (the API space ndim). */
function ndimDimension(d: number): Dimension {
  const t = ND_TEXT[d];
  return {
    id: `${d}-in-${d}`,
    label: `${d}D in ${d}D`,
    summary: `Equal ${d}D polytopes or balls of diameter 1 packed in the smallest copy of a ${d}D polytope or ball.`,
    families: [
      {
        key: `ndim-${d}`,
        space: "ndim",
        scope: { field: "d", at: 1, values: [String(d)] },
        start: `ndim/${d}/cube/cube/12`,
        name: "Polytopes and balls",
        short: "Polytopes",
        summary: `Pieces and container are each one of ${t.shapes}, or a ball: any pair, up to ${t.maxP} pieces.`,
        available: true,
      },
    ],
  };
}

/** The d-cubes' P ranges (docs/spaces/BOUNDARY.md 4B). */
const NCUBE_MAX_P: Record<number, number> = { 5: 160, 6: 384, 7: 400, 8: 400, 9: 400, 10: 400 };

/** "4D on 5D" to "9D on 10D": balls and unit cubes on the boundary of the d-cube, d = 5 to 10 (the API space boundary). */
function ncubeDimension(d: number): Dimension {
  return {
    id: `${d - 1}-on-${d}`,
    label: `${d - 1}D on ${d}D`,
    summary: `Balls of diameter 1 or unit ${d - 1}-cubes on the boundary of the ${d}-cube, with no container. The boundary is a closed ${d - 1}D space made of ${2 * d} cells, each a ${d - 1}-cube, glued face to face.`,
    families: [
      {
        key: `boundary-${d}cube`,
        space: "boundary",
        scope: { field: "polytope", at: 1, values: [`${d}cube`] },
        start: `boundary/${d}cube/ball/9`,
        name: `${d}-cube boundary`,
        short: `${d}-cube`,
        summary: `Balls or axis-aligned unit cubes on the boundary of the ${d}-cube of edge s, up to ${NCUBE_MAX_P[d]} pieces. The score is s.`,
        available: true,
      },
    ],
  };
}

const BASE_DIMENSIONS: readonly Dimension[] = [
  {
    id: "2-in-2",
    label: "2D in 2D",
    summary: "Polygons and disks inside a polygon or disk container, in the flat plane or at constant curvature.",
    families: [
      {
        space: "polygon",
        name: "Flat",
        short: "Flat",
        summary:
          "Polygons of side 1 or disks of diameter 1 in the smallest regular polygon or disk in the flat plane: the board's classic problems.",
        available: true,
      },
      {
        space: "curved",
        name: "Curved",
        short: "Curved",
        summary:
          "Pieces in a container at constant curvature K: in the hyperbolic plane for K < 0, on a sphere for K > 0.",
        available: true,
      },
    ],
  },
  {
    id: "2-on-3",
    label: "2D on 3D",
    summary:
      "Polygons and disks on a closed surface, with no container: the surface of a polyhedron (a Platonic, Archimedean, Catalan or Johnson solid, a prism or an antiprism), a sphere or a flat torus. A flat torus is pictured as a surface in 3D, although it cannot be embedded there both flat and smooth.",
    families: [
      {
        space: "surface",
        name: "Polyhedral surfaces",
        short: "Surfaces",
        summary:
          "Pieces on the surface of one of 141 solids: a Platonic, Archimedean, Catalan or Johnson solid, a prism or an antiprism, with or without folding over its edges.",
        available: true,
      },
      {
        space: "sphere",
        name: "Sphere",
        short: "Sphere",
        summary: "Pieces on a whole sphere: the score is the sphere's radius.",
        available: true,
      },
      {
        space: "torus",
        name: "Torus",
        short: "Torus",
        summary: "Pieces on a flat torus, a plane that repeats along a lattice.",
        available: true,
      },
    ],
  },
  {
    id: "3-in-3",
    label: "3D in 3D",
    summary: "Polyhedra packed inside a polyhedron, balls packed inside a ball or a solid, and solids packed inside a ball.",
    families: [
      {
        space: "polyhedra",
        name: "Polyhedra",
        short: "Polyhedra",
        summary:
          "Equal polyhedra in the smallest copy of a container polyhedron: Platonic, Archimedean, Catalan and Johnson solids, prisms and antiprisms, up to 500 pieces.",
        available: true,
      },
      {
        space: "balls",
        name: "Balls",
        short: "Balls",
        summary:
          "Balls of diameter 1 in the smallest ball or in one of 186 solids, and each of those solids in the smallest ball, up to 100 pieces.",
        available: true,
      },
    ],
  },
  {
    id: "3-on-4",
    label: "3D on 4D",
    summary:
      "Solids on the boundary of a 4D polytope, with no container. The boundary is a closed 3D space made of solid cells, as the surface of a solid is made of flat faces.",
    families: [
      {
        space: "boundary",
        scope: {
          field: "polytope",
          at: 1,
          values: ["tesseract", "5cell", "16cell", "24cell", "120cell", "600cell"],
        },
        start: "boundary/tesseract/ball/9",
        name: "Polytope boundaries",
        short: "Boundaries",
        summary:
          "Balls of diameter 1 or the Platonic solids of edge 1 on the boundary of a 4D polytope: the tesseract, the 4D cube, whose boundary is eight cubes glued face to face, the 5-cell, 16-cell, 24-cell, 120-cell and 600-cell.",
        available: true,
      },
    ],
  },
];

/**
 * Every dimension, in the order of the pieces' dimension and, for the same one, "in" before "on":
 * 2D in 2D, 2D on 3D, 3D in 3D, 3D on 4D, 4D in 4D, 4D on 5D, 5D in 5D, ... 9D on 10D, 10D in 10D.
 */
export const DIMENSIONS: readonly Dimension[] = [
  ...BASE_DIMENSIONS,
  ...[4, 5, 6, 7, 8, 9, 10].flatMap((d) => (d === 10 ? [ndimDimension(d)] : [ndimDimension(d), ncubeDimension(d + 1)])),
];

/** One group of the board's Problem set menu: a dimension, and an option per family. */
export interface MenuGroup {
  label: string;
  options: { text: string; value: string; disabled: boolean }[];
}

/**
 * The board's Problem set menu: a group per dimension, an option per family, in the table's order.
 * A family that is not available yet is listed, disabled, as "coming soon".
 */
export const menuGroups = (dimensions: readonly Dimension[] = DIMENSIONS): MenuGroup[] =>
  dimensions.map((d) => ({
    label: d.label,
    options: d.families.map((f) => ({
      text: `${d.label}: ${f.name}${f.available ? "" : " (coming soon)"}`,
      value: familyKey(f),
      disabled: !f.available,
    })),
  }));

/** The families of a dimension that can be shown now. */
export const availableFamilies = (d: Dimension) =>
  d.families.filter((f): f is Extract<Family, { available: true }> => f.available);

/** A dimension can be chosen when at least one of its families is available. */
export const isAvailable = (d: Dimension): boolean => d.families.some((f) => f.available);

/** What names a family in the menu and in code: its key, else its space. */
export const familyKey = (f: Pick<Family, "key" | "space">): string => f.key ?? f.space;

/** The scope's value in a problem id ("5" in "ndim/5/cube/ball/3"), or null. */
export const scopeValue = (scope: FamilyScope, id: string): string | null => id.split("/")[scope.at] ?? null;

/** True when a scope holds the problem id (its field is one of the scope's values); no scope holds every id of the space. */
export function scopeHasId(scope: FamilyScope | null | undefined, id: string): boolean {
  if (!scope) return true;
  const v = scopeValue(scope, id);
  return v !== null && scope.values.includes(v);
}

/** True when the family shows the problem id. */
export const familyHasId = (f: Family, id: string): boolean => scopeHasId(f.scope, id);

/**
 * The family a board address or menu choice names: the available family with this key; else of this space, the one
 * whose scope holds the id (when an id is given), else the space's first available family. Null when none shows it.
 */
export function whereIs(target: string, id?: string | null): { dimension: Dimension; family: Family } | null {
  return find(target, id, true);
}

/** The same for a family that may be unavailable (a dimension that is coming soon). */
export function familyOfSpace(target: string, id?: string | null) {
  return find(target, id, false);
}

function find(target: string, id: string | null | undefined, onlyAvailable: boolean) {
  const all = DIMENSIONS.flatMap((dimension) =>
    dimension.families.filter((f) => f.available || !onlyAvailable).map((family) => ({ dimension, family })),
  );
  // a problem id names its family by the scope that holds it; else the menu's key; else the space's first family
  return (
    (id ? all.find((x) => x.family.space === target && x.family.scope && familyHasId(x.family, id)) : undefined) ??
    all.find((x) => familyKey(x.family) === target) ??
    all.find((x) => x.family.space === target) ??
    null
  );
}

/**
 * The family of a records table's state: the one whose scope holds the filter on the scope's field; the space's first
 * family when the filter is not a single value (or the space has one family).
 */
export function familyOfFilters(space: string, filters: Record<string, string>) {
  let first: { dimension: Dimension; family: Family } | null = null;
  for (const dimension of DIMENSIONS)
    for (const family of dimension.families) {
      if (family.space !== space) continue;
      first ??= { dimension, family };
      const v = family.scope ? filters[family.scope.field] : undefined;
      if (family.scope && v !== undefined && family.scope.values.includes(v)) return { dimension, family };
    }
  return first;
}

/** The table filters that open a family's records: its scope's field fixed to the one value, none for a scope of several. */
export const familyFilters = (f: Family): Record<string, string> =>
  f.scope && f.scope.values.length === 1 ? { [f.scope.field]: f.scope.values[0] } : {};

/** "2D on 3D, Sphere": a family with its dimension, for sentences. */
export function familyTitle(target: string, id?: string | null): string {
  const where = familyOfSpace(target, id);
  return where ? `${where.dimension.label}, ${where.family.name}` : target;
}
