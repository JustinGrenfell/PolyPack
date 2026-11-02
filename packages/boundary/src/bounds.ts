/**
 * Lower bounds on s for balls of diameter 1 on the boundary of a regular 4-polytope of edge s, and the facts about each
 * container that they rest on (docs/spaces/BOUNDARY.md section 5).
 *
 *  - Volume: developments are local isometries, injective on interiors, with disjoint images, so P balls of volume pi/6
 *    fit only if P pi / 6 <= V s^3, V the volume of the boundary at s = 1 (cells x the cell's volume).
 *  - One piece: a ball's centre needs distance >= 1/2 from every hinge, and every point of a regular cell lies within rho
 *    of an edge of the cell, rho the distance from the cell's centre to its edges (the midradius of the cell). Proof: a
 *    point p of the cell lies on a segment from the centre O to a point q of a facet F of the cell; q lies within the
 *    inradius r_F of a regular polygon F of the boundary of F, i.e. within r_F of an edge e of F; the distance to the
 *    convex set e is convex, so d(p, e) <= (1 - t) d(O, e) + t d(q, e) <= max(rho, r_F) = rho, as rho^2 = h^2 + r_F^2 for
 *    the distance h of O to F's plane. So s rho >= 1/2.
 *
 * The tesseract's bounds are exact (kernel-validate.ts boundaryBound). For the other containers the floats below are
 * rounded down by a relative 1e-15 (far above the rounding of the few operations; the exact check in
 * python/boundary/general.py tests every stored bound against the exact value).
 */
import { NCUBE_MAX_P, ncubeCheckMs } from "./ncube/bounds";
import { NCUBE_NAMES, isNCubeName } from "./ncube/lattice";
import { ndCheckMs } from "./nd/bounds";
import { ND_MAX_P, ND_NAMES, isNDName } from "./nd/names";
import type { ContainerName, PieceName, PolytopeName } from "./types";

const PHI = (1 + Math.sqrt(5)) / 2;

/** What the bounds need of a container, at unit edge. */
export interface ContainerFacts {
  cells: number;
  /** The cell: tetrahedron, cube, octahedron or dodecahedron. */
  cell: "tetrahedron" | "cube" | "octahedron" | "dodecahedron";
  /** Volume of one cell, and of the whole boundary. */
  cellVolume: number;
  boundaryVolume: number;
  /** Distance from the cell's centre to its edges (rho above) and to its facet planes (the inradius). */
  midradius: number;
  inradius: number;
  /** Cells around each edge and the total angle there, in degrees (the full turn minus the deficit). */
  cellsAtEdge: number;
  coneAngle: number;
}

const CELL = {
  tetrahedron: { volume: Math.sqrt(2) / 12, midradius: 1 / (2 * Math.SQRT2), inradius: Math.sqrt(6) / 12 },
  cube: { volume: 1, midradius: Math.SQRT1_2, inradius: 0.5 },
  octahedron: { volume: Math.SQRT2 / 3, midradius: 0.5, inradius: Math.sqrt(6) / 6 },
  dodecahedron: {
    volume: (15 + 7 * Math.sqrt(5)) / 4,
    midradius: (PHI * PHI) / 2,
    inradius: Math.sqrt(250 + 110 * Math.sqrt(5)) / 20,
  },
} as const;

const SHAPE: Record<
  PolytopeName,
  { cell: ContainerFacts["cell"]; cells: number; cellsAtEdge: number; coneAngle: number }
> = {
  "5cell": { cell: "tetrahedron", cells: 5, cellsAtEdge: 3, coneAngle: 211.58633809652784 },
  tesseract: { cell: "cube", cells: 8, cellsAtEdge: 3, coneAngle: 270 },
  "16cell": { cell: "tetrahedron", cells: 16, cellsAtEdge: 4, coneAngle: 282.1151174620371 },
  "24cell": { cell: "octahedron", cells: 24, cellsAtEdge: 3, coneAngle: 328.4136619034721 },
  "120cell": { cell: "dodecahedron", cells: 120, cellsAtEdge: 3, coneAngle: 349.6951535312339 },
  "600cell": { cell: "tetrahedron", cells: 600, cellsAtEdge: 5, coneAngle: 352.6438968275464 },
};

export function containerFacts(name: PolytopeName): ContainerFacts {
  const sh = SHAPE[name];
  const c = CELL[sh.cell];
  return {
    cells: sh.cells,
    cell: sh.cell,
    cellVolume: c.volume,
    boundaryVolume: sh.cells * c.volume,
    midradius: c.midradius,
    inradius: c.inradius,
    cellsAtEdge: sh.cellsAtEdge,
    coneAngle: sh.coneAngle,
  };
}

/** The proven lower bounds on s of P balls on a container other than the tesseract, and the larger. */
export function ballBounds(
  name: PolytopeName,
  P: number,
): { volume: number; onePiece: number; value: number; kind: "volume" | "one piece" } {
  const f = containerFacts(name);
  const down = 1 - 1e-15;
  const volume = Math.cbrt((P * Math.PI) / 6 / f.boundaryVolume) * down;
  const onePiece = (1 / (2 * f.midradius)) * down;
  return volume > onePiece
    ? { volume, onePiece, value: volume, kind: "volume" }
    : { volume, onePiece, value: onePiece, kind: "one piece" };
}

/** Which pieces each container takes (balls on all six; the tesseract's axis-aligned cubes besides). */
export const CONTAINER_PIECES: Record<ContainerName, readonly PieceName[]> = {
  // the d-cubes (d = 5..10) take balls and axis-aligned cubes
  ...(Object.fromEntries(NCUBE_NAMES.map((n) => [n, ["ball", "cube"] as readonly PieceName[]])) as Record<
    (typeof NCUBE_NAMES)[number],
    readonly PieceName[]
  >),
  // the simplex, orthoplex and demicube (d = 5..10) take balls only
  ...(Object.fromEntries(ND_NAMES.map((n) => [n, ["ball"] as readonly PieceName[]])) as Record<(typeof ND_NAMES)[number], readonly PieceName[]>),
  tesseract: ["ball", "cube"],
  "5cell": ["ball"],
  "16cell": ["ball"],
  "24cell": ["ball"],
  "120cell": ["ball"],
  "600cell": ["ball"],
};

/** The largest P of each container in the catalogue: the tesseract's 16, 24 for the others. */
export const CONTAINER_MAX_P: Record<ContainerName, number> = {
  ...NCUBE_MAX_P,
  ...(Object.fromEntries(ND_NAMES.map((n) => [n, ND_MAX_P])) as Record<(typeof ND_NAMES)[number], number>),
  tesseract: 16,
  "5cell": 24,
  "16cell": 24,
  "24cell": 24,
  "120cell": 24,
  "600cell": 24,
};

/**
 * Milliseconds of validation per ball (validateBoundary at tolerance 1e-9, one machine core, [CALC] 2026-11-01): the measured cost of
 * loose random layouts at P = 24 divided by 24 (tesseract 0.3 from BOUNDARY.md section 12), times 3 for tight packings, whose
 * balls see more chains; raised after review to cover 24-ball crowds measured at 20 ms (5-cell) and 40 ms (600-cell). The chain stage dominates, so a container with small cells relative to the ball costs most.
 */
export const BOUNDARY_BALL_MS: Record<PolytopeName, number> = {
  tesseract: 0.3,
  "5cell": 1.0,
  "16cell": 0.8,
  "24cell": 2.3,
  "600cell": 2.0,
  "120cell": 12,
};

/** Estimated milliseconds of the server's check of a configuration of P pieces. Cubes, which have parts instead of chains, cost 0.3 each. */
export function boundaryCheckMs(problem: { polytope: ContainerName; piece: PieceName; P: number }): number {
  if (isNDName(problem.polytope)) return ndCheckMs(problem.polytope, problem.P);
  if (isNCubeName(problem.polytope)) return ncubeCheckMs(problem.polytope, problem.piece as "ball" | "cube", problem.P);
  return problem.piece === "cube" ? 0.3 * problem.P : BOUNDARY_BALL_MS[problem.polytope] * problem.P;
}

/** A check of at least this many milliseconds waits its turn in the server's heavy queue (the 120-cell from P = 5, the 24-cell from P = 22). */
export const BOUNDARY_HEAVY_MS = 50;

export function boundaryIsHeavy(problem: { polytope: ContainerName; piece: PieceName; P: number }): boolean {
  return boundaryCheckMs(problem) >= BOUNDARY_HEAVY_MS;
}
