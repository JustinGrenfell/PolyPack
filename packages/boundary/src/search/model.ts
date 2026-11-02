/**
 * What the boundary search's driver (search.ts) needs from a kind of piece: a penalty energy of a
 * state with its gradient, the largest violation (for the final scaling), starts and hops.
 *
 * A state holds each piece's anchor cell and centre (3 local coordinates at the state's edge s).
 * The driver moves centres through the charts (chart.ts) and asks the model only for energies at
 * actual positions, so a model never deals with re-anchoring.
 *
 * Models: balls (balls.ts, any polytope; exact lattice words on the tesseract) and axis-aligned
 * cubes on the tesseract (cubes.ts). Turned cubes and other polyhedral pieces fit the same
 * interface with a rotation per piece: State gains a rotation (e.g. a quaternion per piece, carried
 * across a ridge by the gluings' linear part M, which the chart already returns), the driver adds
 * three rotation variables per piece, and the energy is a separating-axis penetration depth of the
 * pieces' parts in each cell (parts.ts partPolyhedron, polyhedron.ts satShape/satDepth, which give
 * the depth and the axis; the gradient follows from the axis and the two supporting vertices), or
 * an LP over the parts' half-spaces.
 */
import type { BoundaryConfig, BoundaryPiece } from "../types";

export interface State {
  cells: Int32Array;
  /** 3 numbers per piece: its centre in its anchor cell's local frame, at the state's edge. */
  U: Float64Array;
  /**
   * Models of rotated pieces only (`rotates`): 9 numbers per piece, its orientation, an orthogonal matrix (row-major) in its anchor's local
   * frame (improper after a gluing that reverses handedness; solids.ts).
   */
  R?: Float64Array;
}

/**
 * What the generic methods' kernel (../kernel.ts) asks of a model beyond the search's energy.
 * Without options, terms() gives exactly what energy() and violation() give.
 */
export interface TermOptions {
  /** Per-piece weights (default all 1): a pair term is scaled by (w_i + w_j) / 2, a piece's own terms by w_i. Energy only. */
  w?: Float64Array | null;
  /** Only the terms of this piece: its own and its pairs (default -1: every term). */
  focus?: number;
  /** When given, ds[0] receives dE/ds with the centres U held fixed. Energy only. */
  ds?: Float64Array | null;
}

export interface PieceModel {
  readonly P: number;
  /**
   * True for models of rotated pieces (solids.ts): the state carries orientations, the search adds three rotation variables per piece after
   * the 3P displacements, and energy() / terms() write a gradient of 6P numbers: the centres' 3P, then the 3P rotation increments theta
   * (R' = exp(theta) R, theta in the anchor's frame).
   */
  readonly rotates?: boolean;
  /** A relaxation counts as feasible below this energy. */
  readonly target: number;
  /**
   * Smallest edge worth trying: a proven lower bound for balls; for aligned cubes the family's
   * own floor (no aligned unit cube fits in a cell of edge < 1) or the volume bound.
   */
  readonly floor: number;
  /** Longest move of the centres in one line-search step. */
  readonly maxStep: number;
  /** Penalty energy at edge s (sum of squared violations) and its gradient with respect to U. */
  energy(st: State, s: number, g: Float64Array): number;
  /** Largest violation (a length; <= 0 when the state is valid by the model's own exact rule). */
  violation(st: State, s: number): number;
  /**
   * The energy (g, when given, receives its gradient with respect to U) or, with `violation`, the
   * largest violation, of the terms the options select, weighted as they say.
   */
  terms(st: State, s: number, g: Float64Array | null, violation: boolean, o: TermOptions): number;
  /** Edge of the random starts. */
  startEdge(): number;
  randomState(s: number, rand: () => number): State;
  /** A copy of the state (at edge s) with one to three pieces moved at random. */
  hop(st: State, s: number, rand: () => number): State;
  /** A constructed start (cubes: the two rings), or null. */
  constructedStart(rand: () => number): { st: State; s: number } | null;
  /** The state of a configuration, or null when the model cannot represent it (a turned cube). */
  stateOf(config: BoundaryConfig): State | null;
  /** The pieces of a state (at edge s), clamped into their closed cells. */
  pieces(st: State, s: number): BoundaryPiece[];
}

export function copyState(st: State): State {
  return st.R ? { cells: Int32Array.from(st.cells), U: Float64Array.from(st.U), R: Float64Array.from(st.R) } : { cells: Int32Array.from(st.cells), U: Float64Array.from(st.U) };
}

/** The state scaled by k about the cells' centres (for a change of edge s -> k s). */
export function scaleState(st: State, k: number): State {
  const U = new Float64Array(st.U.length);
  for (let i = 0; i < U.length; i++) U[i] = st.U[i] * k;
  return st.R ? { cells: Int32Array.from(st.cells), U, R: Float64Array.from(st.R) } : { cells: Int32Array.from(st.cells), U };
}
