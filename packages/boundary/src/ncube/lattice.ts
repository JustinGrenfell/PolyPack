/**
 * The boundary of the d-cube [-s/2, s/2]^d, d = 5..10 ("4D on 5D" to "9D on 10D"): lattice words, the closed-form
 * clearance and the geodesic search by reach (docs/spaces/BOUNDARY.md section 4B).
 *
 * Conventions are the tesseract's with 4 replaced by d (section 2): cell c = 2 i + (sigma > 0) is the facet
 * x_i = sigma s/2, a (d-1)-cube with local coordinates u = the other d - 1 coordinates in increasing order, the cell being
 * [-h, h]^n with h = s/2 and n = d - 1. Every ridge of the d-cube is a right angle, so the cells a straight path crosses
 * unfold into the cubes of the lattice 2h Z^n (the anchor at the origin); the hinges, the (d-3)-faces, are the lattice
 * points where at least two coordinates are odd multiples of h (their codimension in R^n is 2). A straight path is
 * described by the lattice planes y_k = (2m + 1) h it crosses in order (a word). Each word has an end cell and an exact
 * map y -> x = A y + h b from the lattice into R^d, A a signed permutation (stored per lattice axis as a row and a sign)
 * and b an integer vector, built by composing the unfoldings: crossing in direction tau e_k from cell (i, sg), where
 * A e_k = g e_j, enters the cell (j, tj), tj = tau g, by the rotation R of the plane (e_i, e_j) about the ridge,
 * R(sg e_i) = tj e_j and R(tj e_j) = -sg e_i; then A <- R A and b <- R (b - q) + q, q = sg e_i + tj e_j.
 *
 * Completeness (proved in BOUNDARY.md 4.1 and 4B): a straight path of length < 2 s crosses at most two planes of each
 * lattice axis (three consecutive planes are 4 h = 2 s apart in that coordinate), each axis in one direction. Balls need
 * only paths shorter than 1 (their diameter), and a ball fits only if s >= 1/sqrt 2, so 2 s > 1 and the words with
 * |m_k| <= 2 contain every path that matters.
 *
 * Reach. The words are not listed (109,601 per cell for d = 5 and 9,864,101 for d = 6 with two crossings per axis, 633 and
 * 6,331 with one; 10^16 for d = 10) but found for one pair of points at a time, by two exact bounds on the length of any
 * word whose straight path from u joins u to the image y of the target:
 *  1. the path crosses, on every axis k with m_k != 0, the planes up to the last one, at distance
 *     delta_k = (2 |m_k| - 1) h - sign(m_k) u_k from u, in a different coordinate for each axis, so
 *     length^2 >= sum_k delta_k^2 (the support search below is pruned by it);
 *  2. the image is y = 2 h m + yhat, where yhat is a signed permutation of the target's local coordinates, so with
 *     c = 2 h m - u, length^2 = sum_r (c_r + yhat_r)^2 >= sum_r (|c|_(r) - |u2|_(r))^2, the two lists sorted (the
 *     rearrangement inequality). It holds for every order of the crossings, so a support (m) that fails it needs no orders.
 * Only the supports that pass both are expanded into orders of crossings (the multinomial many), and each ordered word
 * is tested exactly as the tesseract's: its image in its lattice cube and its planes crossed in order within [0, 1].
 * The search is capped by a work budget (WordLimitError); the cap makes the check sound but incomplete: a configuration
 * that needs more work is refused, never accepted.
 */

export const NCUBE_NAMES = ["5cube", "6cube", "7cube", "8cube", "9cube", "10cube"] as const;
export type NCubeName = (typeof NCUBE_NAMES)[number];

/** The ambient dimension d of a d-cube name ("7cube" -> 7), or 0 when the name is none. */
export function ncubeDim(name: unknown): number {
  const i = (NCUBE_NAMES as readonly unknown[]).indexOf(name);
  return i < 0 ? 0 : i + 5;
}

export function isNCubeName(x: unknown): x is NCubeName {
  return ncubeDim(x) > 0;
}

const AXES = new Map<number, number[]>();
/** The local axes of cell (i, .) in R^d: every axis except i, increasing. */
export function axesOf(d: number, i: number): number[] {
  const key = 16 * d + i;
  let a = AXES.get(key);
  if (!a) {
    a = [];
    for (let j = 0; j < d; j++) if (j !== i) a.push(j);
    AXES.set(key, a);
  }
  return a;
}

/** The R^d point of local coordinates u (read from U at offset o) of cell c, at half-edge h. */
export function embedCell(d: number, c: number, U: ArrayLike<number>, o: number, h: number): Float64Array {
  const x = new Float64Array(d);
  const i = c >> 1;
  x[i] = c & 1 ? h : -h;
  const ax = axesOf(d, i);
  for (let k = 0; k < d - 1; k++) x[ax[k]] = U[o + k];
  return x;
}

/**
 * Geodesic distance from the local point u (read from U at offset o) to the hinges, in closed form: the nearest hinge of the
 * development is an edge-of-cell hinge of the anchor, a pair of axes a < b with distance sqrt((h - |u_a|)^2 + (h - |u_b|)^2),
 * so the least is the root of the two smallest squared slacks. (Proof: BOUNDARY.md 4.2; the odd multiple of h nearest
 * u_a in [-h, h] is +-h.) For n = d - 1 >= 2 axes.
 */
export function ncubeClearance(U: ArrayLike<number>, o: number, n: number, h: number): number {
  let a1 = Infinity;
  let a2 = Infinity;
  for (let k = 0; k < n; k++) {
    const a = h - Math.abs(U[o + k]);
    if (a < a1) {
      a2 = a1;
      a1 = a;
    } else if (a < a2) a2 = a;
  }
  return Math.hypot(a1, a2);
}

// ------------------------------------------------------------------------------------------ the development of a word

/** The development map of a word so far: end cell, and per lattice axis the row and sign of A's column, and b. */
export interface WordMap {
  cell: number;
  rows: Int8Array;
  sgn: Int8Array;
  b: Int8Array;
}

/** The empty word's map at cell c0: lattice axis k is the cell's k-th local axis. */
export function startMap(d: number, c0: number): WordMap {
  const i = c0 >> 1;
  const n = d - 1;
  const rows = new Int8Array(n);
  const sgn = new Int8Array(n).fill(1);
  const ax = axesOf(d, i);
  for (let k = 0; k < n; k++) rows[k] = ax[k];
  const b = new Int8Array(d);
  b[i] = c0 & 1 ? 1 : -1;
  return { cell: c0, rows, sgn, b };
}

/** The map after crossing the plane of lattice axis k in direction tau (+1 or -1). */
export function stepMap(d: number, w: WordMap, k: number, tau: number): WordMap {
  const i = w.cell >> 1;
  const sg = w.cell & 1 ? 1 : -1;
  const j = w.rows[k];
  const tj = tau * w.sgn[k];
  const n = d - 1;
  const rows = Int8Array.from(w.rows);
  const sgn = Int8Array.from(w.sgn);
  // the column that lay along e_j now lies along e_i
  for (let r = 0; r < n; r++)
    if (w.rows[r] === j) {
      rows[r] = i;
      sgn[r] = -sg * tj * w.sgn[r];
    }
  // b <- R (b - q) + q, q = sg e_i + tj e_j
  const b = Int8Array.from(w.b);
  const wi = w.b[i] - sg;
  const wj = w.b[j] - tj;
  b[i] = -sg * tj * wj + sg;
  b[j] = sg * tj * wi + tj;
  return { cell: 2 * j + (tj > 0 ? 1 : 0), rows, sgn, b };
}

/**
 * The map of a word from the target's side: y_r = M u2 + h t, with M a signed permutation, for the target (cell
 * `w.cell`) with local coordinates u2. perm[r] is the index k of u2 that y_r reads, sign[r] its sign, t[r] the shift in
 * units of h. (The columns of A are tangent to the end cell, so none lies along its normal.)
 */
export function targetMap(d: number, w: WordMap, perm: Uint8Array, sign: Int8Array, t: Float64Array): void {
  const ax = axesOf(d, w.cell >> 1);
  const n = d - 1;
  for (let r = 0; r < n; r++) {
    const row = w.rows[r];
    perm[r] = ax.indexOf(row);
    sign[r] = w.sgn[r];
    t[r] = -w.sgn[r] * w.b[row];
  }
}

// ------------------------------------------------------------------------------------------ the search

/** The work of a check is capped; a search that needs more throws this and the caller refuses the configuration. */
export class WordLimitError extends Error {
  constructor(readonly limit: number) {
    super(`the geodesic search needed more than ${limit} steps`);
  }
}

/** A work budget shared by the searches of one check (steps left). */
export interface WordBudget {
  left: number;
  /** The budget's size, for messages. */
  limit: number;
}

export function newBudget(limit: number): WordBudget {
  return { left: limit, limit };
}

/** One realised straight geodesic from a point to (an image of) another, with the map that produced its image. */
export interface Realised {
  /** Its length. */
  len: number;
  /** The image y of the target in the start's development (n numbers). */
  y: Float64Array;
  /** y_r = sign_r u2[perm_r] + h t_r (sign 0: a constant, the coordinate along the target's own normal). */
  perm: Uint8Array;
  sign: Int8Array;
  t: Float64Array;
}

export interface RealisedOptions {
  /** Only paths shorter than this are looked for (default 1, a ball's diameter). */
  reach?: number;
  /** Skip the empty word (for loops, where the target is the start). */
  excludeEmpty?: boolean;
  /** Slack for a point lying on a plane or a face (default 1e-12 max(1, h)). */
  eps?: number;
}

function ascending(a: Float64Array): Float64Array {
  for (let i = 1; i < a.length; i++) {
    const v = a[i];
    let j = i - 1;
    while (j >= 0 && a[j] > v) {
      a[j + 1] = a[j];
      j--;
    }
    a[j + 1] = v;
  }
  return a;
}

/**
 * All realised straight geodesics, shorter than the reach, from the point (c1, U1[o1..]) to the point (c2, U2[o2..]) at
 * edge s, one per realised word that the bounds below do not rule out. A step whose plane holds both ends (two points of
 * one ridge) puts no constraint on the order, as in the tesseract's tables. Caution as there: such a word can also be a path
 * through a hinge, never shorter than the two clearances together, so for balls whose clearances are >= 1/2 it is >= 1 and
 * changes no verdict.
 *
 * The target is the R^d point x2 (the cells that contain it are the end cells). Three exact bounds prune the search, all
 * valid for every word with the support m at hand (the crossings per lattice axis):
 *  1. the support bound: length^2 >= sum_k delta_k^2, delta_k the distance from u to the last plane crossed on axis k, plus,
 *     for an axis k with m_k = 0, the squared difference (x2[row_k] - u_k)^2: a column of A that is never crossed keeps its
 *     row, and b keeps its value there, so y_k = x2[row_k] is known before any crossing;
 *  2. the assignment bound of the header (rearrangement inequality);
 *  3. during the orders: the image coordinate y_k of an axis is known as soon as its last crossing is made (the crossing
 *     changes only column k of A and the entries of b along the old normal and the old row, and no later crossing of another
 *     axis touches them), so delta_k = y_k - u_k replaces the plane distance in the sum, the image must lie in its lattice
 *     cube, and the crossing times t = (plane - u_k)/(y_k - u_k) known so far must lie in [0, 1] and increase along the
 *     order. A prefix that fails any of them has no realised extension.
 */
export function realisedWords(
  d: number,
  s: number,
  c1: number,
  U1: ArrayLike<number>,
  o1: number,
  c2: number,
  U2: ArrayLike<number>,
  o2: number,
  budget: WordBudget,
  opts: RealisedOptions = {},
): Realised[] {
  const n = d - 1;
  const h = s / 2;
  const reach = opts.reach ?? 1;
  const reach2 = reach * reach;
  const eps = opts.eps ?? 1e-12 * Math.max(1, h);
  const excludeEmpty = opts.excludeEmpty === true;
  const u1 = new Float64Array(n);
  for (let k = 0; k < n; k++) u1[k] = U1[o1 + k];
  const x2 = embedCell(d, c2, U2, o2, h);
  // the cells that hold x2: its anchor and, for a coordinate on a ridge, the neighbour across
  const ends = new Set<number>([c2]);
  for (let t = 0; t < d; t++) if (h - Math.abs(x2[t]) <= eps) ends.add(2 * t + (x2[t] > 0 ? 1 : 0));
  const i2 = c2 >> 1;
  const ax2 = axesOf(d, i2);
  const eAbs = ascending(Float64Array.from(ax2, (a) => Math.abs(x2[a])));
  const start = startMap(d, c1);
  // The normal axis of the end cell is the original row of a crossed axis, or the start's own normal axis (a crossing swaps
  // the normal with the row of the axis it crosses, so the normal after any word is one of the rows of the crossed axes, or
  // the original normal when a swap is undone). So unless an end cell has the start's normal axis, some axis whose original
  // row is an end cell's normal axis must be crossed: when only one row is possible, that axis is crossed in every word.
  const i1 = c1 >> 1;
  let mustCross = -1;
  let mayStay = false;
  const needRows = new Set<number>();
  for (const e of ends) {
    if (e >> 1 === i1) mayStay = true;
    else needRows.add(e >> 1);
  }
  if (!mayStay && needRows.size === 1) mustCross = start.rows.indexOf([...needRows][0]);
  // y_k for an axis never crossed: x2 along its original row
  const known0 = new Float64Array(n);
  for (let k = 0; k < n; k++) known0[k] = x2[start.rows[k]] - u1[k];
  // distance from u to the last plane crossed on the axis, for m = 1, 2 on the + and - sides
  const dPos = new Float64Array(2 * n);
  const dNeg = new Float64Array(2 * n);
  for (let k = 0; k < n; k++)
    for (let m = 1; m <= 2; m++) {
      dPos[2 * k + m - 1] = Math.max(0, (2 * m - 1) * h - u1[k]);
      dNeg[2 * k + m - 1] = Math.max(0, (2 * m - 1) * h + u1[k]);
    }
  const out: Realised[] = [];
  const m = new Int8Array(n);
  const tmp = new Float64Array(n);
  const y = new Float64Array(n);
  const stepK = new Int8Array(2 * n);
  const stepLevel = new Int8Array(2 * n);
  const stepT = new Float64Array(2 * n);
  const taken = new Int8Array(n);
  const rem = new Int8Array(n);
  const yKnown = new Uint8Array(n);
  const planeDist2 = new Float64Array(n);

  /** The time at which the segment u -> y crosses a plane, or NaN when the step puts no constraint (both ends on the plane). */
  const crossingTime = (k: number, level: number, yk: number): number => {
    const plane = h * level;
    if (Math.abs(u1[k] - plane) <= eps && Math.abs(yk - plane) <= eps) return NaN;
    return (plane - u1[k]) / (yk - u1[k]);
  };

  /** Tests the finished word map against the end cells. */
  const finish = (w: WordMap, steps: number) => {
    if (!ends.has(w.cell)) return;
    for (let r = 0; r < n; r++) {
      const row = w.rows[r];
      y[r] = w.sgn[r] * (x2[row] - h * w.b[row]);
      if (!(Math.abs(y[r] - 2 * h * m[r]) <= h + eps)) return;
    }
    let prev = -Infinity;
    for (let q = 0; q < steps; q++) {
      const k = stepK[q];
      const plane = h * stepLevel[q];
      if (Math.abs(u1[k] - plane) <= eps && Math.abs(y[k] - plane) <= eps) continue;
      const t = (plane - u1[k]) / (y[k] - u1[k]);
      if (!(t >= -1e-12 && t <= 1 + 1e-12 && t > prev)) return;
      prev = t;
    }
    let l2 = 0;
    for (let r = 0; r < n; r++) l2 += (y[r] - u1[r]) * (y[r] - u1[r]);
    const real: Realised = {
      len: Math.sqrt(l2),
      y: Float64Array.from(y),
      perm: new Uint8Array(n),
      sign: new Int8Array(n),
      t: new Float64Array(n),
    };
    for (let r = 0; r < n; r++) {
      const row = w.rows[r];
      if (row === i2) {
        // the coordinate along the target's own normal: a constant (sigma h)
        real.sign[r] = 0;
        real.t[r] = (w.sgn[r] * (x2[row] - h * w.b[row])) / h;
      } else {
        real.perm[r] = ax2.indexOf(row);
        real.sign[r] = w.sgn[r];
        real.t[r] = -w.sgn[r] * w.b[row];
      }
    }
    out.push(real);
  };

  /** Expands one support into the orders of its crossings; `partial` is the bound 3 so far (squared). */
  const orders = (w: WordMap, steps: number, total: number, partial: number) => {
    if (steps === total) {
      finish(w, steps);
      return;
    }
    for (let k = 0; k < n; k++) {
      if (rem[k] === 0) continue;
      budget.left--;
      if (budget.left < 0) throw new WordLimitError(budget.limit);
      const tau = m[k] > 0 ? 1 : -1;
      const level = tau * (2 * taken[k] + 1);
      stepK[steps] = k;
      stepLevel[steps] = level;
      stepT[steps] = NaN;
      rem[k]--;
      taken[k]++;
      const nw = stepMap(d, w, k, tau);
      let ok = true;
      let np = partial;
      let setT = false;
      if (rem[k] === 0) {
        // the last crossing of the axis: its image coordinate is known
        const row = nw.rows[k];
        const yk = nw.sgn[k] * (x2[row] - h * nw.b[row]);
        if (!(Math.abs(yk - 2 * h * m[k]) <= h + eps)) ok = false;
        else {
          const delta = yk - u1[k];
          np = partial - planeDist2[k] + delta * delta;
          if (!(np < reach2)) ok = false;
          else {
            // times of this axis's crossings (steps so far), then the order among all known times
            for (let q = 0; q <= steps; q++) if (stepK[q] === k) stepT[q] = crossingTime(k, stepLevel[q], yk);
            setT = true;
            let prev = -Infinity;
            for (let q = 0; q <= steps && ok; q++) {
              const t = stepT[q];
              if (Number.isNaN(t)) continue;
              if (!(t >= -1e-12 && t <= 1 + 1e-12 && t >= prev - 1e-12)) ok = false;
              prev = Math.max(prev, t);
            }
          }
        }
      }
      if (ok) orders(nw, steps + 1, total, np);
      if (setT) for (let q = 0; q <= steps; q++) if (stepK[q] === k) stepT[q] = NaN;
      taken[k]--;
      rem[k]++;
    }
  };

  /** Chooses the crossings of axis r onward (the support), pruning by bound 1. */
  const support = (r: number, cost: number, total: number) => {
    budget.left--;
    if (budget.left < 0) throw new WordLimitError(budget.limit);
    if (r === n) {
      if (total === 0 && excludeEmpty) return;
      if (!mayStay) {
        let reaches = false;
        for (let q = 0; q < n && !reaches; q++) if (m[q] !== 0 && needRows.has(start.rows[q])) reaches = true;
        if (!reaches) return;
      }
      // bound 2, the assignment bound
      for (let q = 0; q < n; q++) tmp[q] = Math.abs(2 * h * m[q] - u1[q]);
      ascending(tmp);
      let lb = 0;
      for (let q = 0; q < n; q++) lb += (tmp[q] - eAbs[q]) * (tmp[q] - eAbs[q]);
      if (!(lb < reach2)) return;
      for (let q = 0; q < n; q++) {
        rem[q] = Math.abs(m[q]);
        yKnown[q] = 0;
        planeDist2[q] = m[q] === 0 ? 0 : (m[q] > 0 ? dPos : dNeg)[2 * q + Math.abs(m[q]) - 1] ** 2;
      }
      taken.fill(0);
      orders(start, 0, total, cost);
      return;
    }
    // no crossing on this axis: y_r is known
    const c0 = cost + known0[r] * known0[r];
    if (c0 < reach2 && r !== mustCross) {
      m[r] = 0;
      support(r + 1, c0, total);
    }
    for (let mm = 1; mm <= 2; mm++) {
      const cp = cost + dPos[2 * r + mm - 1] ** 2;
      if (cp < reach2) {
        m[r] = mm;
        support(r + 1, cp, total + mm);
      }
      const cn = cost + dNeg[2 * r + mm - 1] ** 2;
      if (cn < reach2) {
        m[r] = -mm;
        support(r + 1, cn, total + mm);
      }
    }
    m[r] = 0;
  };
  support(0, 0, 0);
  return out;
}
