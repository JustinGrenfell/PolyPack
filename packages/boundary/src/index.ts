/** The boundary space: (d-1)-dimensional pieces on the boundary of a convex d-polytope ("3D on 4D"). */
export * from "./types";
export { BALL_RADIUS } from "./units";
export { ISO_IDENTITY, isoApply, isoCompose, isoFromAffine, isoInverse, isoLinear, isoToAffine, matrixToQuaternion, quaternionToMatrix } from "./linalg";
export type { Iso3, V3 } from "./linalg";
export { convexHull } from "./hull";
export type { HullFacet, HullRidge } from "./hull";
export { boundaryVolume, buildPolytope, cellExcess, compareNormals, facetPlane, fromAmbient, getPolytope, isPolytopeName, tesseractVertices, toAmbient } from "./polytope";
export { BOUNDARY_DOMAIN_MAX_P, BOUNDARY_PROBLEM_COUNT, boundaryDomain } from "./domain";
export { BOUNDARY_BALL_MS, BOUNDARY_HEAVY_MS, CONTAINER_MAX_P, CONTAINER_PIECES, ballBounds, boundaryCheckMs, pieceBounds, boundaryIsHeavy, containerFacts } from "./bounds";
export type { ContainerFacts } from "./bounds";
export { POLYTOPE_PIECES, isPolytopePiece, pieceShape, rotationOnto } from "./pieces";
export type { PieceShape } from "./pieces";
export { RAW_EDGE, rawVertices, unitVertices } from "./containers";
export { MAX_P, formatBoundaryId, isPieceCount, isPieceName, parseBoundaryId, standardPiece } from "./ids";
export {
  ChainLimitError,
  MIN_WINDOW_ANGLE,
  RIDGE_SNAP,
  clampToCell,
  clipSegmentToWindow,
  enumerateChains,
  geodesicDistance,
  geodesicLoop,
  hingeClearance,
  inWindow,
  pathLengths,
  traceGeodesic,
} from "./chains";
export type { BoundaryPoint, ChainNode, ChainOptions } from "./chains";
export { CUBE_FACES, CUBE_VERTICES, PART_CLIP_EPS, PART_EPS, developPiece, partMesh, partPolyhedron, pieceParts, piecePose, pieceReach, windowInCell } from "./parts";
export type { Development } from "./parts";
export { clipPolyhedron, polyhedronFrom, polyhedronThickness, polyhedronVertices, polyhedronVolume, satDepth, satShape } from "./polyhedron";
export type { Polyhedron, SatShape } from "./polyhedron";
export { bfsNet, defaultNet, netFromTree } from "./net";
export { MAX_S, validateBoundary, validateOn } from "./validate";
export type { ValidateOptions } from "./validate";
export { TesseractWords, embed as tesseractEmbed, tesseractClearance } from "./tesseract";
export * from "./search";
export {
  NCUBE_NAMES,
  WordLimitError,
  axesOf as ncubeAxes,
  embedCell as ncubeEmbed,
  isNCubeName,
  ncubeClearance,
  ncubeDim,
  realisedWords as ncubeRealisedWords,
  newBudget as ncubeBudget,
} from "./ncube/lattice";
export { NCUBE_MAX_P, NCUBE_TITLES, ncubeBound, ncubeCheckMs, ncubeProblemBound, unitBallVolume } from "./ncube/bounds";
export { NCUBE_MAX_WORK, cubeParts as ncubeCubeParts, validateNCube } from "./ncube/validate";
export type { NCubeProblem } from "./ncube/validate";
export { NCUBE_SERVER_RULE, ncubeFormErrors, ncubeServerRule } from "./ncube/rule";
export { ncubeKernel } from "./ncube/kernel";
export { searchNCube } from "./ncube/search";
export { cellLayout as ncubeCellLayout, ringLayout as ncubeRingLayout } from "./ncube/seeds";
export { isContainerName } from "./polytope";
export * from "./kernel";
