/**
 * Vertex lists of the six regular convex 4-polytopes, the containers of the boundary space, with unit edge and centred at the
 * origin (the contract's convention, types.ts). Each list is generated in a fixed order that python/boundary/general.py
 * repeats in exact arithmetic (the exact vertices of a container are the ones here, up to the scale), so a vertex index means
 * the same vertex on both sides; a test compares them.
 *
 *  - tesseract: polytope.ts (tesseractVertices).
 *  - 5-cell: the regular simplex (1, 1, 1, -1/sqrt 5), (1, -1, -1, -1/sqrt 5), (-1, 1, -1, -1/sqrt 5), (-1, -1, 1, -1/sqrt 5),
 *    (0, 0, 0, 4/sqrt 5), edge sqrt 8.
 *  - 16-cell: +-e_i, edge sqrt 2.
 *  - 24-cell: the permutations of (+-1, +-1, 0, 0), edge sqrt 2.
 *  - 600-cell: unit circumradius; (+-1, 0, 0, 0) and permutations, (+-1/2, +-1/2, +-1/2, +-1/2), and the even permutations of
 *    (+-phi, +-1, +-1/phi, 0) / 2; edge 1/phi.
 *  - 120-cell: the standard coordinates of circumradius 2 sqrt 2 and edge 3 - sqrt 5 (Coxeter): the permutations of
 *    (0, 0, +-2, +-2), (+-1, +-1, +-1, +-sqrt 5), (+-phi^-2, +-phi, +-phi, +-phi), (+-phi^-1, +-phi^-1, +-phi^-1, +-phi^2), and
 *    the even permutations of (0, +-phi^-2, +-1, +-phi^2), (0, +-phi^-1, +-phi, +-sqrt 5), (+-phi^-1, +-1, +-phi, +-2).
 * The generators below return the vertices before the unit-edge scale (RAW) so that the exact side can repeat them without
 * a square root; `unitVertices` scales.
 */
import type { PolytopeName } from "./types";

const PHI = (1 + Math.sqrt(5)) / 2;
const SQ5 = Math.sqrt(5);

/** The 12 even permutations of (0, 1, 2, 3), in lexicographic order. */
export const EVEN_PERMUTATIONS: readonly (readonly number[])[] = (() => {
  const out: number[][] = [];
  const perm = (a: number[], rest: number[]) => {
    if (!rest.length) {
      let inv = 0;
      for (let i = 0; i < 4; i++) for (let j = i + 1; j < 4; j++) if (a[i] > a[j]) inv++;
      if (inv % 2 === 0) out.push(a);
      return;
    }
    for (const r of rest)
      perm(
        [...a, r],
        rest.filter((x) => x !== r),
      );
  };
  perm([], [0, 1, 2, 3]);
  return out;
})();

/** All 24 permutations of (0, 1, 2, 3), in lexicographic order. */
const ALL_PERMUTATIONS: readonly (readonly number[])[] = (() => {
  const out: number[][] = [];
  const perm = (a: number[], rest: number[]) => {
    if (!rest.length) return void out.push(a);
    for (const r of rest)
      perm(
        [...a, r],
        rest.filter((x) => x !== r),
      );
  };
  perm([], [0, 1, 2, 3]);
  return out;
})();

/**
 * The points obtained from `base` by the given permutations (point[perm[i]] = base[i]) and every sign choice on the nonzero
 * entries, deduplicated in order of first appearance. A sign pattern runs through the nonzero entries of `base` in order,
 * the first entry's sign changing slowest, + before -.
 */
export function signedPoints(base: readonly number[], perms: readonly (readonly number[])[]): number[][] {
  const nz = base.map((x, i) => (x !== 0 ? i : -1)).filter((i) => i >= 0);
  const out: number[][] = [];
  const seen = new Set<string>();
  for (const perm of perms)
    for (let m = 0; m < 1 << nz.length; m++) {
      const p = [0, 0, 0, 0];
      base.forEach((x, i) => {
        const k = nz.indexOf(i);
        const sign = k >= 0 && (m >> (nz.length - 1 - k)) & 1 ? -1 : 1;
        p[perm[i]] = sign * x;
      });
      const key = p.join(",");
      if (!seen.has(key)) {
        seen.add(key);
        out.push(p);
      }
    }
  return out;
}

/** Raw (unscaled) vertices of a container other than the tesseract, in the order the exact side repeats. */
export function rawVertices(name: Exclude<PolytopeName, "tesseract">): number[][] {
  switch (name) {
    case "5cell": {
      const a = 1 / SQ5;
      return [
        [1, 1, 1, -a],
        [1, -1, -1, -a],
        [-1, 1, -1, -a],
        [-1, -1, 1, -a],
        [0, 0, 0, 4 * a],
      ];
    }
    case "16cell": {
      const out: number[][] = [];
      for (let i = 0; i < 4; i++)
        for (const sg of [-1, 1]) {
          const p = [0, 0, 0, 0];
          p[i] = sg;
          out.push(p);
        }
      return out;
    }
    case "24cell":
      return signedPoints([1, 1, 0, 0], ALL_PERMUTATIONS);
    case "600cell":
      return [
        ...signedPoints([1, 0, 0, 0], ALL_PERMUTATIONS),
        ...signedPoints([0.5, 0.5, 0.5, 0.5], [[0, 1, 2, 3]]),
        ...signedPoints([PHI / 2, 0.5, 1 / PHI / 2, 0], EVEN_PERMUTATIONS),
      ];
    case "120cell": {
      const ip = 1 / PHI;
      return [
        ...signedPoints([0, 0, 2, 2], ALL_PERMUTATIONS),
        ...signedPoints([1, 1, 1, SQ5], ALL_PERMUTATIONS),
        ...signedPoints([ip * ip, PHI, PHI, PHI], ALL_PERMUTATIONS),
        ...signedPoints([ip, ip, ip, PHI * PHI], ALL_PERMUTATIONS),
        ...signedPoints([0, ip * ip, 1, PHI * PHI], EVEN_PERMUTATIONS),
        ...signedPoints([0, ip, PHI, SQ5], EVEN_PERMUTATIONS),
        ...signedPoints([ip, 1, PHI, 2], EVEN_PERMUTATIONS),
      ];
    }
  }
}

/** The edge of the raw vertex list (before scaling), as a float. */
export const RAW_EDGE: Record<Exclude<PolytopeName, "tesseract">, number> = {
  "5cell": Math.sqrt(8),
  "16cell": Math.SQRT2,
  "24cell": Math.SQRT2,
  "600cell": 1 / PHI,
  "120cell": 3 - SQ5,
};

/** Unit-edge vertices of a container other than the tesseract. */
export function unitVertices(name: Exclude<PolytopeName, "tesseract">): number[][] {
  const e = RAW_EDGE[name];
  return rawVertices(name).map((v) => v.map((x) => x / e));
}
