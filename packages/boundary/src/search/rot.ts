/**
 * Rotations of R^3 for the pieces' orientations (search/solids.ts, kernel-solids.ts): 3 x 3 matrices stored row-major in 9
 * numbers, an orientation being an orthogonal matrix (improper after a gluing that reverses handedness: the Platonic solids
 * are all mirror-symmetric, so an improper orientation is a proper one of the same piece, shape.mirror).
 */

export type M9 = number[] | Float64Array;

export function mul9(a: M9, ao: number, b: M9, bo: number, out: M9, oo: number): void {
  const a0 = a[ao], a1 = a[ao + 1], a2 = a[ao + 2], a3 = a[ao + 3], a4 = a[ao + 4], a5 = a[ao + 5], a6 = a[ao + 6], a7 = a[ao + 7], a8 = a[ao + 8];
  const b0 = b[bo], b1 = b[bo + 1], b2 = b[bo + 2], b3 = b[bo + 3], b4 = b[bo + 4], b5 = b[bo + 5], b6 = b[bo + 6], b7 = b[bo + 7], b8 = b[bo + 8];
  out[oo] = a0 * b0 + a1 * b3 + a2 * b6;
  out[oo + 1] = a0 * b1 + a1 * b4 + a2 * b7;
  out[oo + 2] = a0 * b2 + a1 * b5 + a2 * b8;
  out[oo + 3] = a3 * b0 + a4 * b3 + a5 * b6;
  out[oo + 4] = a3 * b1 + a4 * b4 + a5 * b7;
  out[oo + 5] = a3 * b2 + a4 * b5 + a5 * b8;
  out[oo + 6] = a6 * b0 + a7 * b3 + a8 * b6;
  out[oo + 7] = a6 * b1 + a7 * b4 + a8 * b7;
  out[oo + 8] = a6 * b2 + a7 * b5 + a8 * b8;
}

export function det9(m: M9, o = 0): number {
  return m[o] * (m[o + 4] * m[o + 8] - m[o + 5] * m[o + 7]) - m[o + 1] * (m[o + 3] * m[o + 8] - m[o + 5] * m[o + 6]) + m[o + 2] * (m[o + 3] * m[o + 7] - m[o + 4] * m[o + 6]);
}

export const IDENTITY9 = [1, 0, 0, 0, 1, 0, 0, 0, 1];

/** exp of the cross-product matrix of w (Rodrigues). */
export function expRot(w0: number, w1: number, w2: number, out: M9, oo = 0): void {
  const t2 = w0 * w0 + w1 * w1 + w2 * w2;
  const t = Math.sqrt(t2);
  let a: number;
  let b: number;
  if (t < 1e-4) {
    a = 1 - t2 / 6 + (t2 * t2) / 120;
    b = 0.5 - t2 / 24 + (t2 * t2) / 720;
  } else {
    a = Math.sin(t) / t;
    b = (1 - Math.cos(t)) / t2;
  }
  // R = I + a K + b K^2, K = [w]x
  out[oo] = 1 + b * (-w1 * w1 - w2 * w2);
  out[oo + 1] = -a * w2 + b * w0 * w1;
  out[oo + 2] = a * w1 + b * w0 * w2;
  out[oo + 3] = a * w2 + b * w0 * w1;
  out[oo + 4] = 1 + b * (-w0 * w0 - w2 * w2);
  out[oo + 5] = -a * w0 + b * w1 * w2;
  out[oo + 6] = -a * w1 + b * w0 * w2;
  out[oo + 7] = a * w0 + b * w1 * w2;
  out[oo + 8] = 1 + b * (-w0 * w0 - w1 * w1);
}

/** The left Jacobian of the exponential: exp(w + dw) = exp(Jl dw) exp(w) to first order. */
export function leftJacobian(w0: number, w1: number, w2: number, out: M9, oo = 0): void {
  const t2 = w0 * w0 + w1 * w1 + w2 * w2;
  const t = Math.sqrt(t2);
  let b: number;
  let c: number;
  if (t < 1e-3) {
    b = 0.5 - t2 / 24;
    c = 1 / 6 - t2 / 120;
  } else {
    b = (1 - Math.cos(t)) / t2;
    c = (t - Math.sin(t)) / (t2 * t);
  }
  // I + b K + c K^2
  out[oo] = 1 + c * (-w1 * w1 - w2 * w2);
  out[oo + 1] = -b * w2 + c * w0 * w1;
  out[oo + 2] = b * w1 + c * w0 * w2;
  out[oo + 3] = b * w2 + c * w0 * w1;
  out[oo + 4] = 1 + c * (-w0 * w0 - w2 * w2);
  out[oo + 5] = -b * w0 + c * w1 * w2;
  out[oo + 6] = -b * w1 + c * w0 * w2;
  out[oo + 7] = b * w0 + c * w1 * w2;
  out[oo + 8] = 1 + c * (-w0 * w0 - w1 * w1);
}

/** The rotation vector of a proper rotation matrix (angle in [0, pi]). */
export function logRot(m: M9, o = 0): [number, number, number] {
  const c = Math.max(-1, Math.min(1, (m[o] + m[o + 4] + m[o + 8] - 1) / 2));
  const x = m[o + 7] - m[o + 5];
  const y = m[o + 2] - m[o + 6];
  const z = m[o + 3] - m[o + 1];
  const sn = Math.hypot(x, y, z) / 2;
  const ang = Math.atan2(sn, c);
  if (sn < 1e-12) {
    if (c > 0) return [x / 2, y / 2, z / 2];
    // angle pi: the axis from the symmetric part
    const a0 = Math.sqrt(Math.max(0, (m[o] + 1) / 2));
    const a1 = Math.sqrt(Math.max(0, (m[o + 4] + 1) / 2));
    const a2 = Math.sqrt(Math.max(0, (m[o + 8] + 1) / 2));
    const ax: [number, number, number] = [a0, m[o + 1] >= 0 ? a1 : -a1, m[o + 2] >= 0 ? a2 : -a2];
    if (a0 < 1e-9) ax[0] = 0;
    return [ax[0] * Math.PI, ax[1] * Math.PI, ax[2] * Math.PI];
  }
  const k = ang / (2 * sn);
  return [x * k, y * k, z * k];
}

/** Gram-Schmidt on the rows of an almost orthogonal matrix, in place (keeps the handedness). */
export function orthonormalize(m: M9, o = 0): void {
  let n0 = Math.hypot(m[o], m[o + 1], m[o + 2]);
  m[o] /= n0; m[o + 1] /= n0; m[o + 2] /= n0;
  let d = m[o] * m[o + 3] + m[o + 1] * m[o + 4] + m[o + 2] * m[o + 5];
  m[o + 3] -= d * m[o]; m[o + 4] -= d * m[o + 1]; m[o + 5] -= d * m[o + 2];
  n0 = Math.hypot(m[o + 3], m[o + 4], m[o + 5]);
  m[o + 3] /= n0; m[o + 4] /= n0; m[o + 5] /= n0;
  // the third row from the cross product, with the sign of the original third row
  const cx = m[o + 1] * m[o + 5] - m[o + 2] * m[o + 4];
  const cy = m[o + 2] * m[o + 3] - m[o] * m[o + 5];
  const cz = m[o] * m[o + 4] - m[o + 1] * m[o + 3];
  d = cx * m[o + 6] + cy * m[o + 7] + cz * m[o + 8];
  const sg = d < 0 ? -1 : 1;
  m[o + 6] = sg * cx; m[o + 7] = sg * cy; m[o + 8] = sg * cz;
}

/** A uniformly random rotation matrix (random unit quaternion). */
export function randomRotation(rand: () => number): number[] {
  let w: number, x: number, y: number, z: number, n: number;
  do {
    w = 2 * rand() - 1; x = 2 * rand() - 1; y = 2 * rand() - 1; z = 2 * rand() - 1;
    n = w * w + x * x + y * y + z * z;
  } while (n > 1 || n < 1e-6);
  const r = 1 / Math.sqrt(n);
  w *= r; x *= r; y *= r; z *= r;
  return [
    1 - 2 * (y * y + z * z), 2 * (x * y - w * z), 2 * (x * z + w * y),
    2 * (x * y + w * z), 1 - 2 * (x * x + z * z), 2 * (y * z - w * x),
    2 * (x * z - w * y), 2 * (y * z + w * x), 1 - 2 * (x * x + y * y),
  ];
}
