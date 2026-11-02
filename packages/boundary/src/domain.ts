/**
 * The domain of the boundary space for the spaces listing (packages/spaces domain.ts): balls on all six regular 4-polytopes
 * (P = 1 to 24 on the five added 2026-11-01, 1 to 16 on the tesseract) and the tesseract's axis-aligned cubes (P = 1 to 16); since
 * 2026-11-02 the five Platonic solids, freely rotated, on all six (the cube on the tesseract being the aligned one), same ranges; and the d-cubes.
 * The tesseract's 32 problems have research records; the others have none, only their lower bounds.
 */
import { CONTAINER_MAX_P, CONTAINER_PIECES, ballBounds, pieceBounds } from "./bounds";
import { boundaryBound } from "./kernel-validate";
import { NCUBE_TITLES, ncubeProblemBound } from "./ncube/bounds";
import { isNCubeName } from "./ncube/lattice";
import { ndBoundKind, ndProblemBound } from "./nd/bounds";
import { ND_TITLES, isNDName, ndParse } from "./nd/names";
import type { ContainerName, PieceName, PolytopeName, PolytopePieceName } from "./types";
import { CONTAINERS, PIECES, POLYTOPES } from "./types";

const TITLES: Record<ContainerName, string> = {
  ...NCUBE_TITLES,
  ...ND_TITLES,
  tesseract: "Tesseract",
  "5cell": "5-cell",
  "16cell": "16-cell",
  "24cell": "24-cell",
  "120cell": "120-cell",
  "600cell": "600-cell",
};

/** The most pieces of any container (the P field's overall range). */
export const BOUNDARY_DOMAIN_MAX_P = Math.max(...CONTAINERS.map((n) => CONTAINER_MAX_P[n]));

/** Number of problems of the domain. */
export const BOUNDARY_PROBLEM_COUNT = CONTAINERS.reduce(
  (n, p) => n + CONTAINER_PIECES[p].length * CONTAINER_MAX_P[p],
  0,
);

const isPolytope = (x: unknown): x is ContainerName =>
  typeof x === "string" && (CONTAINERS as readonly string[]).includes(x);
const isPiece = (x: unknown): x is PieceName => typeof x === "string" && (PIECES as readonly string[]).includes(x);

/** The bound of a problem of the domain as the listing shows it, or undefined when the id is outside the domain. */
function boundEntry(polytope: ContainerName, piece: PieceName, P: number): { value: number; label: string } | null {
  if (isNCubeName(polytope)) {
    const value = ncubeProblemBound(polytope, piece as "ball" | "cube", P);
    return piece === "ball"
      ? {
          value,
          label:
            "[CALC] larger of the volume bound (P balls of diameter 1 have at most the (d-1)-volume 2d s^(d-1) of the boundary) and one ball (1/sqrt 2: the centre is 1/2 from every (d-3)-face, and no point is farther than s/sqrt 2)",
        }
      : {
          value,
          label:
            "[CALC] larger of the volume bound (P <= 2d s^(d-1)) and one cube (1: with s < 1 it crosses every axis of its cell)",
        };
  }
  if (isNDName(polytope)) {
    const q = ndParse(polytope)!;
    return {
      value: ndProblemBound(polytope, P),
      label:
        ndBoundKind(q.family, q.d, P) === "volume"
          ? "[CALC] volume: the P balls of radius 1/2 are Euclidean balls (the centre is 1/2 from the hinges) with disjoint interiors, so P V_(d-1) 2^-(d-1) <= s^(d-1) times the boundary's (d-1)-volume at unit edge"
          : "[CALC] one ball: its centre must be 1/2 from every (d-3)-face, and every point of a cell lies within rho of a (d-3)-face of the cell, rho the largest distance from the cell's centre to one",
    };
  }
  // the Platonic solids, turned (the tesseract's cube is the aligned one of its records, below): the volume and one-piece bounds of 4C.4
  if (piece !== "ball" && !(piece === "cube" && polytope === "tesseract")) {
    const b = pieceBounds(polytope as PolytopeName, piece as PolytopePieceName, P);
    return {
      value: b.value,
      label:
        b.kind === "volume"
          ? `[CALC] volume: the developed ${piece}s have disjoint interiors, so P V <= V_boundary s^3, V the piece's volume`
          : `[CALC] one piece: it contains its inscribed ball, whose centre must be at distance >= the inradius from the edges, and every point of a cell lies within s rho of the cell's edges, rho the distance from a cell's centre to them`,
    };
  }
  if (polytope === "tesseract") {
    const value = boundaryBound(piece, P);
    return piece === "ball"
      ? {
          value,
          label:
            "[CALC] larger of the volume bound (P pi/6 <= 8 s^3) and one ball (1/sqrt 2: the centre is 1/2 from every edge)",
        }
      : { value, label: "[CALC] larger of the volume bound (P <= 8 s^3) and one cube (1/sqrt 2)" };
  }
  const b = ballBounds(polytope as PolytopeName, P);
  return {
    value: b.value,
    label:
      b.kind === "volume"
        ? "[CALC] volume: the developed balls have disjoint interiors, so P pi/6 <= V s^3, V the volume of the boundary at unit edge"
        : "[CALC] one ball: its centre must be at distance >= 1/2 from the edges, and every point of a cell lies within s rho of the cell's edges, rho the distance from a cell's centre to them",
  };
}

export const boundaryDomain = {
  info: {
    id: "boundary/<polytope>/<piece>/<P>",
    fields: [
      { key: "polytope", label: "Polytope", values: [...CONTAINERS], titles: TITLES },
      {
        key: "piece",
        label: "Piece",
        values: [...PIECES],
        titles: { ball: "Balls", cube: "Cubes", tetrahedron: "Tetrahedra", octahedron: "Octahedra", dodecahedron: "Dodecahedra", icosahedron: "Icosahedra" },
        // a field whose values depend on the polytope: balls and the Platonic solids on the 4-polytopes (the tesseract's cube is the aligned one, the others' a rotated one), balls and aligned cubes on the d-cubes
        dependent: { on: "polytope", values: Object.fromEntries(CONTAINERS.map((n) => [n, [...CONTAINER_PIECES[n]]])) },
      },
      {
        key: "P",
        label: "Pieces (P)",
        min: 1,
        max: BOUNDARY_DOMAIN_MAX_P,
        dependent: {
          on: "polytope",
          ranges: Object.fromEntries(CONTAINERS.map((n) => [n, [1, CONTAINER_MAX_P[n]] as [number, number]])),
        },
      },
    ],
    problems: BOUNDARY_PROBLEM_COUNT,
  },
  bound(id: { fields: Readonly<Record<string, string | number>>; P: number }) {
    const { polytope, piece } = id.fields;
    if (!isPolytope(polytope) || !isPiece(piece) || !CONTAINER_PIECES[polytope].includes(piece)) return undefined;
    if (!Number.isInteger(id.P) || id.P < 1 || id.P > CONTAINER_MAX_P[polytope]) return undefined;
    return boundEntry(polytope, piece, id.P);
  },
  /** Small P first: P = 1, 2, ... over the containers and pieces (the volunteers' removed work list drew in this order). */
  sample(limit: number): string[] {
    const out: string[] = [];
    for (let P = 1; P <= BOUNDARY_DOMAIN_MAX_P; P++)
      for (const n of CONTAINERS)
        for (const piece of CONTAINER_PIECES[n]) {
          if (P > CONTAINER_MAX_P[n]) continue;
          if (out.length >= limit) return out;
          out.push(`boundary/${n}/${piece}/${P}`);
        }
    return out;
  },
};
