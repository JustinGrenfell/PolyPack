/**
 * Moving a centre on the boundary of the d-cube: the exponential map in the cells' local charts (the d-cube twin of
 * search/chart.ts, which is written for 3-dimensional cells).
 *
 * A centre (cell c, local u) moved by a vector v of c's local frame follows the segment u + t v of c's development:
 * straight inside the cube [-h, h]^n; where it reaches a face u_k = +-h it crosses into the neighbouring cell
 * (j, tau), j the R^d axis of local axis k, tau the sign of v_k, unfolded by the rotation about the ridge (lattice.ts):
 * the old normal axis i becomes a local axis of the new cell with coordinate sg h on the ridge, and the direction's
 * component v_k turns into -sg tau v_k along it. Every step is a signed permutation of the local axes, so the linear
 * part of the whole map is one: v_end[r] = sgn[r] v_start[perm[r]]. A segment that reaches a hinge (two coordinates at
 * +-h together) has no continuation: the boundary is not flat there.
 */
import { axesOf } from "./lattice";

/** Crossings after which a walk gives up (a displacement of this many faces is far beyond any geodesic that matters). */
const MAX_CROSSINGS = 64;

/**
 * Walks from (cell c0, start[so..so+n)) by v[vo..vo+n) at edge s. Writes the end point into outU[ou..), the signed
 * permutation of the displacement into perm[po..po+n) and sgn[po..po+n) (when given), and returns the end cell, or -1
 * when the segment meets a hinge or crosses more than MAX_CROSSINGS faces.
 */
export function walkN(
  d: number,
  c0: number,
  start: ArrayLike<number>,
  so: number,
  v: ArrayLike<number>,
  vo: number,
  s: number,
  outU: Float64Array,
  ou: number,
  perm: Int32Array | null,
  sgn: Float64Array | null,
  po: number,
): number {
  const n = d - 1;
  const h = s / 2;
  let p = new Float64Array(n);
  let e = new Float64Array(n);
  let len2 = 0;
  for (let k = 0; k < n; k++) {
    p[k] = start[so + k];
    e[k] = v[vo + k];
    len2 += e[k] * e[k];
  }
  const tp = new Int32Array(n);
  const ts = new Float64Array(n);
  for (let k = 0; k < n; k++) {
    tp[k] = k;
    ts[k] = 1;
  }
  const len = Math.sqrt(len2);
  let c = c0;
  if (!Number.isFinite(len)) return -1;
  if (len > 0) {
    for (let k = 0; k < n; k++) e[k] /= len;
    let rem = len;
    const tolH = 1e-12 * Math.max(1, s);
    const q = new Float64Array(n);
    for (let step = 0; ; step++) {
      if (step > MAX_CROSSINGS) return -1;
      let best = Infinity;
      let kb = -1;
      for (let k = 0; k < n; k++) {
        if (e[k] === 0) continue;
        const t = ((e[k] > 0 ? h : -h) - p[k]) / e[k];
        if (t < best) {
          best = t;
          kb = k;
        }
      }
      if (kb < 0 || best >= rem) {
        for (let k = 0; k < n; k++) p[k] += rem * e[k];
        break;
      }
      if (best < 0) best = 0;
      for (let k = 0; k < n; k++) q[k] = p[k] + best * e[k];
      // the crossing point lies on another face too: a point of a hinge
      for (let k = 0; k < n; k++) if (k !== kb && h - Math.abs(q[k]) < tolH) return -1;
      const i = c >> 1;
      const sg = c & 1 ? 1 : -1;
      const ax = axesOf(d, i);
      const j = ax[kb];
      const tau = e[kb] > 0 ? 1 : -1;
      const ax2 = axesOf(d, j);
      const np = new Float64Array(n);
      const ne = new Float64Array(n);
      const ntp = new Int32Array(n);
      const nts = new Float64Array(n);
      for (let m = 0; m < n; m++) {
        const a = ax2[m];
        if (a === i) {
          np[m] = sg * h;
          ne[m] = -sg * tau * e[kb];
          ntp[m] = tp[kb];
          nts[m] = -sg * tau * ts[kb];
        } else {
          const l = ax.indexOf(a);
          np[m] = q[l];
          ne[m] = e[l];
          ntp[m] = tp[l];
          nts[m] = ts[l];
        }
      }
      p = np;
      e = ne;
      tp.set(ntp);
      ts.set(nts);
      c = 2 * j + (tau > 0 ? 1 : 0);
      rem -= best;
    }
  }
  // rounding may leave the point a hair outside its cell
  for (let k = 0; k < n; k++) outU[ou + k] = Math.max(-h, Math.min(h, p[k]));
  if (perm && sgn) {
    for (let k = 0; k < n; k++) {
      perm[po + k] = tp[k];
      sgn[po + k] = ts[k];
    }
  }
  return c;
}
