/**
 * The simplex, orthoplex and demicube boundaries for d = 5..10 as cell complexes (docs/spaces/BOUNDARY.md section 4D.1).
 *
 * The container is a convex d-polytope P of unit edge centred at its vertex centroid, P = { x : n_c . x <= h_c for every cell c },
 * n_c the unit outward normal of the cell (facet) c. The boundary is the union of the cells; each cell is a (d-1)-dimensional
 * polytope in the hyperplane H_c = { n_c . x = h_c }, glued to its neighbours along ridges. The coordinates are generated in
 * BASE UNITS (edge sqrt 2, rational or, for the demicube's mixed ridges, in Q(sqrt d)) and divided by E0 = sqrt 2:
 *  - simplex (D = d + 1, in the hyperplane sum x = 1): vertices e_0..e_d. Cell i is the facet x_i = 0 (vertices e_j, j != i);
 *    normal (1 - (d+1) e_i) / sqrt(d (d+1)), offset 1 / sqrt(d (d+1)). Every two cells are adjacent. d + 1 cells.
 *  - orthoplex (D = d): vertices +-e_i. Cell c is the facet sigma . x = 1, sigma_i = -1 exactly when bit i of c is set
 *    (vertices sigma_i e_i); normal sigma / sqrt d, offset 1 / sqrt d. Cells c and c xor 2^k are adjacent. 2^d cells.
 *  - demicube (D = d): vertices (+-1/2)^d with an even number of minus signs. Cells 0..2d-1: the facet x_i = a/2 (a = +-1) has
 *    id 2 i + (a > 0), a (d-1)-demicube, normal a e_i, offset 1/2. Cells 2d + k: the facet w . x = (d-2)/2 at an odd-parity sign
 *    vector w (bit i of the mask set exactly when w_i = -1; k the rank of the mask among the masks of odd popcount), a
 *    (d-1)-simplex, normal w / sqrt d, offset (d-2) / (2 sqrt d). 2d + 2^(d-1) cells. A demicube facet x_i = a/2 is adjacent to the
 *    other demicube facets x_j = b/2 (j != i) and to the simplex facets with w_i = a; a simplex facet to its d demicube facets.
 *    [EST: Coxeter, Regular Polytopes, 8.7 for the demicube's facets; CALC: the adjacencies, checked against a convex hull in the tests]
 *
 * Ridges are the (d-2)-faces shared by two adjacent cells; the singular set of the boundary is the (d-3)-faces, the hinges. A
 * hinge of a cell is the intersection of two of its ridges that is (d-3)-dimensional (`hinges`: pairs of indices into `nbr`).
 *
 * Local coordinates of a cell: u in R^n (n = d - 1) at edge s means x = s o_c + sum_k u_k f_k, o_c the cell's centre (the centroid of
 * its vertices) and f_0..f_(n-1) the cell's frame: Gram-Schmidt (run twice) on the projections of e_0, e_1, ... onto the cell's
 * tangent space, skipping a vector whose residual is under 1e-9, until n are found. For the demicube's x_i = a/2 cells this is the
 * other coordinates in increasing order.
 *
 * Everything is built lazily per cell (the 10-orthoplex has 1,024 cells and the 10-demicube 532).
 */
import { type NDFamily, type NDName, ndCellCount, ndParse } from "./names";

/** The edge of the base polytopes before scaling to unit edge. */
export const E0 = Math.SQRT2;

export interface NDCell {
  c: number;
  /** 0: simplex or orthoplex cell; 1: a demicube's (d-1)-demicube facet; 2: a demicube's simplex facet. */
  kind: 0 | 1 | 2;
  /** Unit outward normal, and offset at unit edge (n . x = h on the cell, <= h in the polytope). */
  n: Float64Array;
  h: number;
  /** Neighbouring cells, their unit normals (flat, k * D), offsets, and a point of each ridge at unit edge (flat). */
  nbr: Int32Array;
  nn: Float64Array;
  hn: Float64Array;
  rp: Float64Array;
  /** |nu_k|^2 = 1 - (n_k . n)^2, the squared length of the tangential part of neighbour k's normal. */
  nu2: Float64Array;
  /** The cell's centre at unit edge, and its frame (n orthonormal vectors of R^D). */
  o: Float64Array;
  frame: Float64Array[];
  /** Pairs (k, l) of indices into nbr whose ridges meet in a hinge, flat. */
  hinges: Int32Array;
}

const popcount = (x: number): number => {
  let c = 0;
  for (; x; x &= x - 1) c++;
  return c;
};

export class NDContainer {
  readonly name: NDName;
  readonly family: NDFamily;
  readonly d: number;
  readonly n: number;
  /** Dimension of the ambient space the coordinates live in (d + 1 for the simplex, d otherwise). */
  readonly D: number;
  readonly C: number;
  /** The simplex lies in the hyperplane flat.v . x = flat.b (at unit edge); null for the others. */
  readonly flat: { v: Float64Array; b: number } | null;
  private readonly cells: (NDCell | undefined)[];
  private readonly oddMasks: Int32Array;
  private readonly oddIndex: Int32Array;
  private readonly sd: number;
  private cum: Float64Array | null = null;

  constructor(name: NDName) {
    const p = ndParse(name);
    if (!p) throw new Error(`NDContainer: unknown container ${String(name)}`);
    this.name = name;
    this.family = p.family;
    this.d = p.d;
    this.n = p.d - 1;
    this.D = p.family === "simplex" ? p.d + 1 : p.d;
    this.C = ndCellCount(p.family, p.d);
    this.sd = Math.sqrt(p.d);
    this.cells = new Array(this.C);
    if (p.family === "simplex") {
      const v = new Float64Array(this.D).fill(1 / Math.sqrt(this.D));
      this.flat = { v, b: 1 / (Math.sqrt(this.D) * E0) };
    } else this.flat = null;
    const odd: number[] = [];
    for (let m = 0; m < 2 ** p.d; m++) if (popcount(m) % 2 === 1) odd.push(m);
    this.oddMasks = Int32Array.from(odd);
    this.oddIndex = new Int32Array(2 ** p.d).fill(-1);
    odd.forEach((m, k) => (this.oddIndex[m] = k));
  }

  cell(c: number): NDCell {
    let cl = this.cells[c];
    if (!cl) {
      if (!Number.isInteger(c) || c < 0 || c >= this.C) throw new Error(`cell ${c} is out of range`);
      cl = this.build(c);
      this.cells[c] = cl;
    }
    return cl;
  }

  /** The sign vector w of a demicube simplex cell (id >= 2d), as +-1 entries. */
  private signs(mask: number): Float64Array {
    const w = new Float64Array(this.d);
    for (let i = 0; i < this.d; i++) w[i] = (mask >> i) & 1 ? -1 : 1;
    return w;
  }

  /** Unit normal and unit-edge offset of cell c, and its kind. */
  private normalOf(c: number): { n: Float64Array; h: number; kind: 0 | 1 | 2 } {
    const { d, D, sd } = this;
    const n = new Float64Array(D);
    if (this.family === "simplex") {
      const k = Math.sqrt(d * (d + 1));
      for (let j = 0; j < D; j++) n[j] = 1 / k;
      n[c] = (1 - (d + 1)) / k;
      return { n, h: 1 / k / E0, kind: 0 };
    }
    if (this.family === "orthoplex") {
      for (let i = 0; i < d; i++) n[i] = ((c >> i) & 1 ? -1 : 1) / sd;
      return { n, h: 1 / sd / E0, kind: 0 };
    }
    if (c < 2 * d) {
      n[c >> 1] = c & 1 ? 1 : -1;
      return { n, h: 0.5 / E0, kind: 1 };
    }
    const w = this.signs(this.oddMasks[c - 2 * d]);
    for (let i = 0; i < d; i++) n[i] = w[i] / sd;
    return { n, h: (d - 2) / (2 * sd) / E0, kind: 2 };
  }

  private neighbours(c: number): number[] {
    const { d, family } = this;
    const out: number[] = [];
    if (family === "simplex") {
      for (let g = 0; g <= d; g++) if (g !== c) out.push(g);
    } else if (family === "orthoplex") {
      for (let k = 0; k < d; k++) out.push(c ^ (1 << k));
    } else if (c < 2 * d) {
      const i = c >> 1;
      const neg = !(c & 1);
      for (let j = 0; j < d; j++) if (j !== i) out.push(2 * j, 2 * j + 1);
      for (let k = 0; k < this.oddMasks.length; k++) {
        const m = this.oddMasks[k];
        if (((m >> i) & 1) === (neg ? 1 : 0)) out.push(2 * d + k);
      }
    } else {
      const m = this.oddMasks[c - 2 * d];
      for (let i = 0; i < d; i++) out.push(2 * i + ((m >> i) & 1 ? 0 : 1));
    }
    return out;
  }

  /** A point of the ridge between adjacent cells c and g, at unit edge (the centroid of the ridge's vertices). */
  ridgePoint(c: number, g: number): Float64Array {
    const { d, D, family } = this;
    const x = new Float64Array(D);
    if (family === "simplex") {
      for (let j = 0; j < D; j++) x[j] = 1 / (d - 1);
      x[c] = 0;
      x[g] = 0;
    } else if (family === "orthoplex") {
      const k = Math.log2(c ^ g);
      for (let i = 0; i < d; i++) x[i] = i === k ? 0 : ((c >> i) & 1 ? -1 : 1) / (d - 1);
    } else if (c < 2 * d && g < 2 * d) {
      x[c >> 1] = c & 1 ? 0.5 : -0.5;
      x[g >> 1] = g & 1 ? 0.5 : -0.5;
    } else {
      // a demicube facet x_i = +-1/2 and a simplex facet w: the centroid of the simplex's vertices off coordinate i
      const dem = c < 2 * d ? c : g;
      const sim = c < 2 * d ? g : c;
      const i = dem >> 1;
      const w = this.signs(this.oddMasks[sim - 2 * d]);
      for (let k = 0; k < d; k++) x[k] = w[k] / 2 - (k === i ? 0 : w[k] / (d - 1));
    }
    for (let j = 0; j < D; j++) x[j] /= E0;
    return x;
  }

  private centre(c: number): Float64Array {
    const { d, D, family } = this;
    const o = new Float64Array(D);
    if (family === "simplex") {
      for (let j = 0; j < D; j++) o[j] = 1 / d;
      o[c] = 0;
    } else if (family === "orthoplex") {
      for (let i = 0; i < d; i++) o[i] = ((c >> i) & 1 ? -1 : 1) / d;
    } else if (c < 2 * d) o[c >> 1] = c & 1 ? 0.5 : -0.5;
    else {
      const w = this.signs(this.oddMasks[c - 2 * d]);
      for (let i = 0; i < d; i++) o[i] = (w[i] * (d - 2)) / (2 * d);
    }
    for (let j = 0; j < D; j++) o[j] /= E0;
    return o;
  }

  private build(c: number): NDCell {
    const { D, n: nn_, d } = this;
    const { n, h, kind } = this.normalOf(c);
    const nb = this.neighbours(c);
    const F = nb.length;
    const nn = new Float64Array(F * D);
    const hn = new Float64Array(F);
    const rp = new Float64Array(F * D);
    const nu2 = new Float64Array(F);
    for (let k = 0; k < F; k++) {
      const g = this.normalOf(nb[k]);
      nn.set(g.n, k * D);
      hn[k] = g.h;
      rp.set(this.ridgePoint(c, nb[k]), k * D);
      let dot = 0;
      for (let j = 0; j < D; j++) dot += g.n[j] * n[j];
      nu2[k] = 1 - dot * dot;
    }
    // hinge pairs
    const pairs: number[] = [];
    {
      if (kind === 1) {
        // demicube facet: (dem j, dem j') with j != j' for d >= 6; (dem j, simplex w) with w_j = the sign of the demicube facet
        const nDem = 2 * (d - 1);
        for (let k = 0; k < nDem; k++) {
          for (let l = k + 1; l < nDem; l++) if (d >= 6 && nb[k] >> 1 !== nb[l] >> 1) pairs.push(k, l);
          for (let l = nDem; l < F; l++) {
            const m = this.oddMasks[nb[l] - 2 * d];
            const j = nb[k] >> 1;
            const wj = (m >> j) & 1 ? -1 : 1;
            if ((nb[k] & 1 ? 1 : -1) === wj) pairs.push(k, l);
          }
        }
      } else {
        for (let k = 0; k < F; k++) for (let l = k + 1; l < F; l++) pairs.push(k, l);
      }
    }
    // frame
    const cons: Float64Array[] = [n];
    if (this.flat) cons.push(this.flat.v);
    const frame: Float64Array[] = [];
    const proj = (v: Float64Array, against: Float64Array[]) => {
      for (const b of against) {
        let dt = 0;
        for (let j = 0; j < D; j++) dt += b[j] * v[j];
        for (let j = 0; j < D; j++) v[j] -= dt * b[j];
      }
    };
    for (let k = 0; k < D && frame.length < nn_; k++) {
      const v = new Float64Array(D);
      v[k] = 1;
      for (let pass = 0; pass < 2; pass++) {
        proj(v, cons);
        proj(v, frame);
      }
      let nm = 0;
      for (let j = 0; j < D; j++) nm += v[j] * v[j];
      nm = Math.sqrt(nm);
      if (nm < 1e-9) continue;
      for (let j = 0; j < D; j++) v[j] /= nm;
      frame.push(v);
    }
    if (frame.length !== nn_) throw new Error(`frame of cell ${c} of the ${this.name} has ${frame.length} vectors`);
    return {
      c,
      kind,
      n,
      h,
      nbr: Int32Array.from(nb),
      nn,
      hn,
      rp,
      nu2,
      o: this.centre(c),
      frame,
      hinges: Int32Array.from(pairs),
    };
  }

  /**
   * A random cell, chosen with probability proportional to its (d-1)-volume (the demicube's simplex facets are tiny beside its
   * (d-1)-demicube facets); `rand` returns numbers in [0, 1).
   */
  pickCell(rand: () => number): number {
    if (this.family !== "demicube") return Math.min(this.C - 1, Math.floor(rand() * this.C));
    if (!this.cum) {
      const n = this.n;
      let nf = 1;
      for (let i = 2; i <= n; i++) nf *= i;
      const wSim = Math.sqrt(this.d) / nf;
      const wDem = 1 - 2 ** (n - 1) / nf;
      const w: number[] = [];
      for (let c = 0; c < this.C; c++) w.push(c < 2 * this.d ? wDem : wSim);
      const total = w.reduce((a, b) => a + b, 0);
      let acc = 0;
      this.cum = Float64Array.from(w, (x) => (acc += x / total));
    }
    const r = rand();
    let lo = 0;
    let hi = this.C - 1;
    while (lo < hi) {
      const mid = (lo + hi) >> 1;
      if (this.cum[mid] >= r) hi = mid;
      else lo = mid + 1;
    }
    return lo;
  }

  /** The cell's name for readers: "x3 = 0" (simplex), the sign vector "+-+-+" of the orthoplex facet sigma . x = 1, "+x2" or "-x2" and "w = -++-+" for the demicube's facets. */
  facetLabel(c: number): string {
    const { d, family } = this;
    const sgn = (m: number) => Array.from({ length: d }, (_, i) => ((m >> i) & 1 ? "-" : "+")).join("");
    if (family === "simplex") return `x${c + 1} = 0`;
    if (family === "orthoplex") return sgn(c);
    if (c < 2 * d) return `${c & 1 ? "+" : "-"}x${(c >> 1) + 1}`;
    return `w = ${sgn(this.oddMasks[c - 2 * d])}`;
  }

  /**
   * The container's vertices at edge s in the engine's ambient coordinates: e_j, +-e_i or the even sign vectors / 2, all scaled by s / sqrt 2. The
   * simplex's centroid is (1/D) 1, not the origin: the ambient coordinates are the engine's, centres from toAmbient are in the same ones.
   */
  vertices(s = 1): Float64Array[] {
    const { d, D, family } = this;
    const f = s / E0;
    const out: Float64Array[] = [];
    if (family === "simplex") {
      for (let j = 0; j < D; j++) {
        const v = new Float64Array(D);
        v[j] = f;
        out.push(v);
      }
    } else if (family === "orthoplex") {
      for (let i = 0; i < d; i++)
        for (const sg of [-1, 1]) {
          const v = new Float64Array(D);
          v[i] = sg * f;
          out.push(v);
        }
    } else {
      for (let m = 0; m < 2 ** d; m++) {
        if (popcount(m) % 2) continue;
        const v = new Float64Array(D);
        for (let i = 0; i < d; i++) v[i] = ((m >> i) & 1 ? -0.5 : 0.5) * f;
        out.push(v);
      }
    }
    return out;
  }

  /** The ambient point of local coordinates u (n numbers, at edge s) of cell c. */
  toAmbient(
    c: number,
    u: ArrayLike<number>,
    s: number,
    uo = 0,
    out: Float64Array = new Float64Array(this.D),
  ): Float64Array {
    const cl = this.cell(c);
    const { D, n } = this;
    for (let j = 0; j < D; j++) out[j] = s * cl.o[j];
    for (let k = 0; k < n; k++) {
      const uk = u[uo + k];
      const f = cl.frame[k];
      for (let j = 0; j < D; j++) out[j] += uk * f[j];
    }
    return out;
  }

  /** The local coordinates of an ambient point of cell c's hyperplane (the inverse of toAmbient). */
  toLocal(c: number, x: ArrayLike<number>, s: number, out: Float64Array = new Float64Array(this.n)): Float64Array {
    const cl = this.cell(c);
    const { D, n } = this;
    for (let k = 0; k < n; k++) {
      const f = cl.frame[k];
      let dt = 0;
      for (let j = 0; j < D; j++) dt += f[j] * (x[j] - s * cl.o[j]);
      out[k] = dt;
    }
    return out;
  }
}

const REGISTRY = new Map<string, NDContainer>();

/** The container of a name (built once). */
export function getND(name: NDName): NDContainer {
  let c = REGISTRY.get(name);
  if (!c) {
    c = new NDContainer(name);
    REGISTRY.set(name, c);
  }
  return c;
}
