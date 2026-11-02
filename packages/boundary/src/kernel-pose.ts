/**
 * Random poses of the boundary kernel (kernel.ts) as displacements. A method undoes a rejected
 * move by writing the piece's x back, so randomPose leaves the anchors alone: it draws a target
 * (a cell, and a point of it) and writes the displacement from the piece's anchor along the
 * shortest geodesic the lattice words give. On the tesseract the geodesics from a point are its lattice
 * words (tesseract.ts): each realised word gives an image y of the target in the anchor's
 * development, the segment from the anchor to y is a straight path avoiding the hinges, and the
 * chart's move along y - anchor ends on the target. Any of them encodes the target exactly, so a
 * uniform target gives a uniform pose; the shortest keeps the displacement short. [CALC] For
 * 80,000 random pairs of points at s = 1 to 3 a word was always found, and the move along it
 * always ended on the target.
 *
 * Any other polytope: the same with the chains of chains.ts instead of words. The target is developed through the chains
 * from the anchor at growing reach (from the chord in R^4, a lower bound of the geodesic, up to a few edges); each chain
 * ending in the target's cell whose closed window holds the target's image gives a straight path avoiding the hinges
 * (the lemma of chains.ts), and the shortest is taken. The enumeration is capped (POSE_CHAINS), so a far target on a large
 * container with small cells gives up (null), and the kernel then draws the end of a random geodesic instead (kernel.ts randomPose).
 */
import { ChainLimitError, clampToCell, enumerateChains, inWindow } from "./chains";
import { V3, isoApply } from "./linalg";
import { TesseractWords } from "./tesseract";
import type { Polytope } from "./types";

let WORDS: TesseractWords | null = null;

/** The lattice words, built once per thread (about 16,000). */
function words(): TesseractWords {
  WORDS ??= new TesseractWords(2);
  return WORDS;
}

/**
 * The displacement (in the frame of `cell`) along the shortest realised geodesic from the point a
 * of `cell` to the point u of `target`, at edge s; null when no word with at most two crossings
 * per lattice axis is realised.
 */
export function displacementTo(
  cell: number,
  a: number[],
  target: number,
  u: number[],
  s: number,
  poly?: Polytope,
): [number, number, number] | null {
  if (poly && poly.name !== "tesseract") return displacementOnChains(poly, cell, a, target, u, s);
  let best: [number, number, number] | null = null;
  let bestLength = Infinity;
  for (const y of words().images(cell, a, target, u, s)) {
    const v: [number, number, number] = [y[0] - a[0], y[1] - a[1], y[2] - a[2]];
    const length = Math.hypot(v[0], v[1], v[2]);
    if (length < bestLength) {
      best = v;
      bestLength = length;
    }
  }
  return best;
}

/** Most chains one pose enumerates (a far target of a container with small cells gives up). */
export const POSE_CHAINS = 1200;

/** Chains version of displacementTo for any polytope. */
function displacementOnChains(
  poly: Polytope,
  cell: number,
  a: number[],
  target: number,
  u: number[],
  s: number,
): [number, number, number] | null {
  const from = clampToCell(poly, cell, a, s);
  const to = clampToCell(poly, target, u, s);
  const xa = ambient(poly, cell, from, s);
  const xb = ambient(poly, target, to, s);
  const chord = Math.hypot(xa[0] - xb[0], xa[1] - xb[1], xa[2] - xb[2], xa[3] - xb[3]);
  // the geodesic is at least the chord; a few edges reach every point of the boundary the catalogue has
  for (let reach = Math.max(1.02 * chord, 0.25 * s); reach <= 6 * s; reach *= 1.6) {
    let nodes;
    try {
      nodes = enumerateChains(poly, s, cell, from, reach, { maxChains: POSE_CHAINS });
    } catch (e) {
      if (e instanceof ChainLimitError) return null;
      throw e;
    }
    let best: [number, number, number] | null = null;
    let bestLength = Infinity;
    for (const node of nodes) {
      if (node.cell !== target) continue;
      const y = isoApply(node.toAnchor, to);
      if (!inWindow(node, from, y)) continue;
      const v: V3 = [y[0] - from[0], y[1] - from[1], y[2] - from[2]];
      const length = Math.hypot(v[0], v[1], v[2]);
      if (length < bestLength) {
        best = v;
        bestLength = length;
      }
    }
    // a path no longer than the reach is the shortest one found; a longer one may have a shorter rival beyond it
    if (best && bestLength <= reach) return best;
  }
  return null;
}

function ambient(poly: Polytope, cell: number, u: number[], s: number): number[] {
  const c = poly.cells[cell];
  const x = c.origin.map((o) => s * o);
  c.frame.forEach((f, k) => f.forEach((fi, i) => (x[i] += u[k] * fi)));
  return x;
}
