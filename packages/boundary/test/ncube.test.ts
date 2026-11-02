/**
 * The boundary of the d-cube, d = 5..10 ("4D on 5D" to "9D on 10D"): the lattice words and the reach search, the closed-form
 * clearance, the float validator (balls and axis-aligned cubes), the bounds, the kernel contract, the seeds and the server rule.
 * The exact certifier (python/boundary/ncube.py) has its own tests and reads the engine's verdicts written by
 * scripts/gen-ncube-crosscheck.ts; here the engine reads the certifier's test packings (python/boundary/fixtures/ncube_packings.json).
 */
import { readFileSync } from "node:fs";
import { fileURLToPath } from "node:url";
import { checkKernel } from "@polypack/methods/contract";
import { mulberry32, squeeze } from "@polypack/methods";
import { describe, expect, it } from "vitest";
import {
  NCUBE_MAX_P,
  NCUBE_MAX_WORK,
  NCUBE_NAMES,
  boundaryBound,
  boundaryDomain,
  boundaryIsHeavy,
  boundaryKernel,
  ncubeAxes,
  ncubeBound,
  ncubeBudget,
  ncubeCellLayout,
  ncubeClearance,
  ncubeDim,
  ncubeRealisedWords,
  ncubeRingLayout,
  ncubeServerRule,
  parseBoundaryId,
  validateBoundary,
  validateNCube,
} from "../src/index";
import { startMap, stepMap } from "../src/ncube/lattice";
import type { WordMap } from "../src/ncube/lattice";

type Piece = "ball" | "cube";
interface Cfg {
  s: number;
  pieces: { cell: number; u: number[] }[];
}
const fixtures = (name: string) =>
  JSON.parse(
    readFileSync(fileURLToPath(new URL(`../../../python/boundary/fixtures/${name}`, import.meta.url)), "utf8"),
  );
const packings: { id: string; polytope: string; piece: Piece; P: number; kind: string; config: Cfg }[] =
  fixtures("ncube_packings.json").packings;
const problem = (d: number, piece: Piece, P: number) => ({
  polytope: `${d}cube` as (typeof NCUBE_NAMES)[number],
  piece,
  P,
});
const verdict = (d: number, piece: Piece, c: Cfg, tol = 0) =>
  validateNCube(problem(d, piece, c.pieces.length), c, { tol });
const zeros = (n: number) => new Array<number>(n).fill(0);

/** Visits every word from c0 with at most mm crossings per axis, in one direction per axis. */
function visitWords(
  d: number,
  c0: number,
  mm: number,
  visit: (w: WordMap, m: number[], steps: [number, number][]) => void,
): void {
  const n = d - 1;
  const rec = (w: WordMap, m: number[], steps: [number, number][]) => {
    visit(w, m, steps);
    for (let k = 0; k < n; k++)
      for (const tau of [-1, 1]) {
        if (m[k] !== 0 && m[k] > 0 !== tau > 0) continue;
        if (Math.abs(m[k] + tau) > mm) continue;
        const m2 = m.slice();
        m2[k] += tau;
        rec(stepMap(d, w, k, tau), m2, [...steps, [k, 2 * m[k] + tau]]);
      }
  };
  rec(startMap(d, c0), zeros(n), []);
}

/** The words, stored (for the smaller enumerations). */
function words(d: number, c0: number, mm: number): { w: WordMap; m: number[]; steps: [number, number][] }[] {
  const out: { w: WordMap; m: number[]; steps: [number, number][] }[] = [];
  visitWords(d, c0, mm, (w, m, steps) => out.push({ w, m, steps }));
  return out;
}

const countFormula = (n: number, mm: number) => {
  const fact = (k: number): number => (k <= 1 ? 1 : k * fact(k - 1));
  const comb = (a: number, b: number) => fact(a) / (fact(b) * fact(a - b));
  if (mm === 1)
    return Array.from({ length: n + 1 }, (_, k) => comb(n, k) * 2 ** k * fact(k)).reduce((a, b) => a + b, 0);
  return Array.from({ length: 2 * n + 1 }, (_, k) => fact(2 * n) / fact(2 * n - k)).reduce((a, b) => a + b, 0);
};

describe("d-cube lattice words", () => {
  it("count 109,601 per cell for d = 5 with two crossings per axis and 633 and 6,331 with one (the research note's)", () => {
    expect([3, 4, 5].map((n) => countFormula(n, 2))).toEqual([1957, 109601, 9864101]);
    expect([3, 4, 5].map((n) => countFormula(n, 1))).toEqual([79, 633, 6331]);
    expect(words(5, 0, 2)).toHaveLength(109601);
    expect(words(5, 6, 1)).toHaveLength(633);
    expect(words(6, 3, 1)).toHaveLength(6331);
  });

  it("develop each cube onto its end cell, and a word that ends in its start cell crosses two axes twice (also for d = 6, all 9,864,101 words)", () => {
    for (const [d, c0] of [
      [5, 0],
      [5, 9],
      [6, 4],
    ] as const) {
      const n = d - 1;
      let own = 0;
      let count = 0;
      visitWords(d, c0, 2, (w, m, steps) => {
        count++;
        const i = w.cell >> 1;
        const sg = w.cell & 1 ? 1 : -1;
        const centre = Array.from(w.b);
        for (let r = 0; r < n; r++) centre[w.rows[r]] += w.sgn[r] * 2 * m[r];
        // the cube of the word lands on its end cell, and the columns of A leave the end cell's normal axis alone
        if (!(centre[i] === sg && centre.every((x, t) => t === i || x === 0)) || w.rows.includes(i))
          throw new Error("not an unfolding");
        if (steps.length && w.cell === c0) {
          own++;
          if (m.filter((x) => Math.abs(x) === 2).length < 2)
            throw new Error("a loop word with fewer than two double crossings");
        }
      });
      expect(count).toBe(countFormula(n, 2));
      expect(own).toBeGreaterThan(0);
    }
  }, 60_000);
});

/** Brute force: the realised words of the full enumeration, as the tesseract's tables test them. */
function bruteLengths(
  d: number,
  s: number,
  c1: number,
  u1: number[],
  c2: number,
  u2: number[],
  mm: number,
  all: ReturnType<typeof words>,
  loop: boolean,
): number[] {
  const n = d - 1;
  const h = s / 2;
  const x2 = new Float64Array(d);
  x2[c2 >> 1] = c2 & 1 ? h : -h;
  ncubeAxes(d, c2 >> 1).forEach((a, k) => (x2[a] = u2[k]));
  const ends = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19].filter(
    (c) => c < 2 * d && Math.abs(x2[c >> 1] - (c & 1 ? h : -h)) <= 1e-12,
  );
  const out: number[] = [];
  for (const { w, m, steps } of all) {
    if (!ends.includes(w.cell) || (loop && steps.length === 0)) continue;
    const y = Array.from({ length: n }, (_, r) => w.sgn[r] * (x2[w.rows[r]] - h * w.b[w.rows[r]]));
    if (!y.every((v, r) => Math.abs(v - 2 * h * m[r]) <= h + 1e-12)) continue;
    let prev = -Infinity;
    let ok = true;
    for (const [k, level] of steps) {
      const plane = h * level;
      if (Math.abs(u1[k] - plane) <= 1e-12 && Math.abs(y[k] - plane) <= 1e-12) continue;
      const t = (plane - u1[k]) / (y[k] - u1[k]);
      if (!(t >= -1e-12 && t <= 1 + 1e-12 && t > prev)) {
        ok = false;
        break;
      }
      prev = t;
    }
    if (ok) out.push(Math.hypot(...y.map((v, r) => v - u1[r])));
  }
  void mm;
  return out.sort((a, b) => a - b);
}

describe("the reach search", () => {
  it("finds exactly the realised words of the full enumeration within the reach (d = 5 with two crossings per axis, d = 6 and 7 with one)", () => {
    const rand = mulberry32(5);
    for (const [d, mm, lo, trials] of [
      [5, 2, 0.75, 60],
      [6, 1, 1, 30],
      [7, 1, 1, 15],
    ] as const) {
      const n = d - 1;
      const fullByCell = new Map<number, ReturnType<typeof words>>();
      let found = 0;
      for (let t = 0; t < trials; t++) {
        let s = 0;
        let c1 = 0;
        let c2 = 0;
        let u1: number[] = [];
        let u2: number[] = [];
        for (;;) {
          s = lo + rand() * (1.6 - lo);
          const h = s / 2;
          c1 = Math.floor(rand() * 2 * d);
          c2 = Math.floor(rand() * 2 * d);
          u1 = Array.from({ length: n }, () => (2 * rand() - 1) * h);
          u2 = Array.from({ length: n }, () => (2 * rand() - 1) * h);
          const U = Float64Array.from([...u1, ...u2]);
          // the chords of R^d: keep pairs that can interact
          const x = (c: number, u: number[]) => {
            const v = new Float64Array(d);
            v[c >> 1] = c & 1 ? h : -h;
            ncubeAxes(d, c >> 1).forEach((a, k) => (v[a] = u[k]));
            return v;
          };
          const a = x(c1, u1);
          const b = x(c2, u2);
          if (Math.hypot(...Array.from(a, (v, i) => v - b[i])) < 1 || U.length === 0) break;
        }
        if (!fullByCell.has(c1)) fullByCell.set(c1, words(d, c1, mm));
        for (const loop of [false, true]) {
          const U = Float64Array.from([...u1, ...(loop ? u1 : u2)]);
          const got = ncubeRealisedWords(d, s, c1, U, 0, loop ? c1 : c2, U, n, ncubeBudget(1e8), { excludeEmpty: loop })
            .map((r) => r.len)
            .filter((l) => l < 1);
          const want = bruteLengths(d, s, c1, u1, loop ? c1 : c2, loop ? u1 : u2, mm, fullByCell.get(c1)!, loop).filter(
            (l) => l < 1,
          );
          expect(got.sort((a, b) => a - b).map((x) => x.toFixed(9))).toEqual(want.map((x) => x.toFixed(9)));
          found += want.length;
        }
      }
      expect(found).toBeGreaterThan(0);
    }
  }, 120_000);

  it("has the closed-form clearance of the nearest hinge lattice point", () => {
    const rand = mulberry32(8);
    for (const d of [5, 7, 10]) {
      const n = d - 1;
      for (let t = 0; t < 60; t++) {
        const h = 0.35 + rand() * 0.9;
        const u = Array.from({ length: n }, () => (2 * rand() - 1) * h);
        let best = Infinity;
        for (let a = 0; a < n; a++)
          for (let b = a + 1; b < n; b++)
            for (let oa = -5; oa <= 5; oa += 2)
              for (let ob = -5; ob <= 5; ob += 2) best = Math.min(best, Math.hypot(u[a] - oa * h, u[b] - ob * h));
        expect(ncubeClearance(u, 0, n, h)).toBeCloseTo(best, 12);
      }
    }
  });

  it("stops at its work budget with an error, never an answer", () => {
    const d = 9;
    const n = d - 1;
    const a = 0.36;
    const U = Float64Array.from([...new Array(n).fill(0.5 - a), ...new Array(n).fill(0.5 - a)]);
    expect(() => ncubeRealisedWords(d, 1, 0, U, 0, 1, U, n, ncubeBudget(3))).toThrow(/needed more than 3 steps/);
  });
});

describe("the d-cube validator", () => {
  it("accepts the certifier's test packings at tolerance 0, refuses their shrunk copies, and states a plain reason", () => {
    expect(packings).toHaveLength(42);
    for (const p of packings) {
      const d = ncubeDim(p.polytope);
      expect(verdict(d, p.piece, p.config).reasons, p.id + " " + p.kind).toEqual([]);
      expect(validateBoundary({ polytope: p.polytope as never, piece: p.piece, P: p.P }, p.config, { tol: 0 }).ok).toBe(
        true,
      );
      const f = 1 - 1e-6;
      const shrunk = {
        s: p.config.s * f,
        pieces: p.config.pieces.map((q) => ({ cell: q.cell, u: q.u.map((x) => x * f) })),
      };
      const v = verdict(d, p.piece, shrunk);
      expect(v.ok, p.id + " " + p.kind).toBe(false);
      expect(v.reasons[0]).toMatch(/^(Ball|Balls|Cube|Cubes) \d|have a total volume/);
    }
  });

  it("accepts the exact tilings: P = 2d k^(d-1) cubes at s = k (d = 5, 6, 8, k = 1, 2) and 160 and 384 balls on the grids at s = 2", () => {
    for (const [d, k] of [
      [5, 1],
      [5, 2],
      [6, 2],
      [8, 1],
    ] as const) {
      const n = d - 1;
      const pts: number[][] = [];
      for (let idx = 0; idx < k ** n; idx++) {
        let r = idx;
        const p: number[] = [];
        for (let a = 0; a < n; a++) {
          p.push((2 * (r % k) - (k - 1)) / 2);
          r = Math.floor(r / k);
        }
        pts.push(p);
      }
      const pieces = Array.from({ length: 2 * d }, (_, c) => pts.map((u) => ({ cell: c, u }))).flat();
      expect(verdict(d, "cube", { s: k, pieces }).reasons).toEqual([]);
      // one cube more, or the tiling a hair smaller, overlaps
      expect(verdict(d, "cube", { s: k, pieces: [...pieces, { cell: 0, u: zeros(n) }] }).ok).toBe(false);
      expect(verdict(d, "cube", { s: k * (1 - 1e-6), pieces }).ok).toBe(false);
    }
    for (const d of [5, 6]) {
      const n = d - 1;
      const pts = Array.from({ length: 2 ** n }, (_, idx) =>
        Array.from({ length: n }, (_, a) => (((idx >> a) & 1) * 2 - 1) / 2),
      );
      const pieces = Array.from({ length: 2 * d }, (_, c) => pts.map((u) => ({ cell: c, u }))).flat();
      expect(verdict(d, "ball", { s: 2, pieces }).reasons).toEqual([]);
      expect(verdict(d, "ball", { s: 2 * (1 - 1e-6), pieces }).ok).toBe(false);
    }
  });

  it("refuses a ball over a (d-3)-face, and keeps the one at 1/2", () => {
    for (let d = 5; d <= 10; d++) {
      const n = d - 1;
      const h = 1.5;
      const ball = (a: number) => ({ s: 3, pieces: [{ cell: 0, u: [h - a, h - a, ...zeros(n - 2)] }] });
      expect(verdict(d, "ball", ball(0.36)).ok).toBe(true);
      const bad = verdict(d, "ball", ball(0.35));
      expect(bad.ok).toBe(false);
      expect(bad.reasons[0]).toMatch(/\(d-3\)-face/);
      expect(verdict(d, "ball", ball(0)).ok).toBe(false);
      // a tolerance forgives 1e-9 of it only
      const edge = Math.sqrt(0.125) - 5e-10;
      expect(verdict(d, "ball", ball(edge), 1e-9).ok).toBe(true);
      expect(verdict(d, "ball", ball(edge), 0).ok).toBe(false);
    }
  });

  it("refuses balls 1e-9 too close across a ridge, with the tolerance deciding the borderline", () => {
    for (let d = 5; d <= 10; d++) {
      const n = d - 1;
      const pair = (eps: number): Cfg => {
        const a = zeros(n);
        a[0] = 1.5 - 0.5;
        const b = zeros(n);
        b[0] = 1.5 - 0.5 + eps;
        return {
          s: 3,
          pieces: [
            { cell: 1, u: a },
            { cell: 3, u: b },
          ],
        };
      };
      // eps > 0 moves B toward the ridge: the unfolded distance is 1 - eps
      expect(verdict(d, "ball", pair(0)).ok).toBe(true);
      expect(verdict(d, "ball", pair(-2e-9)).ok).toBe(true);
      expect(verdict(d, "ball", pair(2e-9)).ok).toBe(false);
      expect(verdict(d, "ball", pair(5e-10), 1e-9).ok).toBe(true);
      expect(verdict(d, "ball", pair(5e-10), 0).ok).toBe(false);
      expect(verdict(d, "ball", pair(2e-9), 1e-9).reasons[0]).toMatch(/^Balls 0 and 1 overlap/);
    }
  });

  it("sees a point on a ridge from both of its cells: the same point listed in the two cells is an overlap, 0.999 along the ridge too, 1 along it is clear", () => {
    for (let d = 5; d <= 10; d++) {
      const n = d - 1;
      const h = 1.5;
      // cell 1 (x0 = h) and cell 3 (x1 = h) share the ridge x0 = x1 = h; the centres lie on it
      const at = (along: number): Cfg => {
        const a = zeros(n);
        a[0] = h;
        const b = zeros(n);
        b[0] = h;
        b[1] = along;
        return {
          s: 3,
          pieces: [
            { cell: 1, u: a },
            { cell: 3, u: b },
          ],
        };
      };
      expect(verdict(d, "ball", at(0)).reasons[0]).toMatch(/their centres are 0 apart/);
      expect(verdict(d, "ball", at(0.999)).ok).toBe(false);
      expect(verdict(d, "ball", at(1)).ok).toBe(true);
    }
  });

  it("finds two balls overlapping around a hinge through a word that crosses two ridges (clearances .6 each, distance .919, the direct way 1.195)", () => {
    // Round the hinge x0 = -h, x1 = h, x2 = h the three cells 0, (1,+), (2,+) fill 270 degrees. A in cell 0 at radius .6 from the hinge and angle
    // 85 degrees (nearly on the ridge to (1,+)), B in cell (2,+) at radius .6 and angle 185 degrees as seen through (1,+): the straight path
    // through (1,+) crosses two ridges and is 2 (.6) sin 50 degrees long; the way across the ridge of cells 0 and (2,+) is 170 degrees round,
    // 2 (.6) sin 85 degrees = 1.195 long. Only the two-crossing word sees the overlap.
    const rad = (x: number) => (x * Math.PI) / 180;
    for (let d = 5; d <= 10; d++) {
      const n = d - 1;
      const h = 1.5;
      const r = 0.6;
      const w = stepMap(d, stepMap(d, startMap(d, 0), 0, 1), 1, 1);
      const u1 = zeros(n);
      u1[0] = h - r * Math.cos(rad(85));
      u1[1] = h - r * Math.sin(rad(85));
      const target = (phi: number) => {
        const y = zeros(n);
        y[0] = h + r * Math.cos(rad(phi - 180));
        y[1] = h + r * Math.sin(rad(phi - 180));
        const x = Array.from(w.b, (b) => h * b);
        for (let k = 0; k < n; k++) x[w.rows[k]] += w.sgn[k] * y[k];
        return ncubeAxes(d, w.cell >> 1).map((a) => x[a]);
      };
      const bad: Cfg = {
        s: 3,
        pieces: [
          { cell: 0, u: u1 },
          { cell: w.cell, u: target(185) },
        ],
      };
      const v = verdict(d, "ball", bad);
      expect(v.ok).toBe(false);
      expect(v.reasons).toHaveLength(1);
      expect(v.reasons[0]).toMatch(/Balls 0 and 1 overlap: their centres are 0\.919\d* apart/);
      // a first-crossing-only validator would call it fine: the other word is 1.195 long
      const fine: Cfg = {
        s: 3,
        pieces: [
          { cell: 0, u: u1 },
          { cell: w.cell, u: target(235) },
        ],
      };
      expect(verdict(d, "ball", fine).reasons).toEqual([]);
    }
  });

  it("refuses a cube over a hinge and cubes overlapping across a ridge", () => {
    for (let d = 5; d <= 10; d++) {
      const n = d - 1;
      const h = 1.5;
      const over = zeros(n);
      over[0] = h - 0.25;
      over[1] = h - 0.25;
      const hinge = verdict(d, "cube", { s: 3, pieces: [{ cell: 0, u: over }] });
      expect(hinge.ok).toBe(false);
      expect(hinge.reasons[0]).toMatch(/\(d-3\)-face/);
      over[1] = 0;
      expect(verdict(d, "cube", { s: 3, pieces: [{ cell: 0, u: over }] }).ok).toBe(true);
      const a = zeros(n);
      a[0] = h - 0.25;
      const pair = (delta: number): Cfg => {
        const b = zeros(n);
        b[0] = h - 0.75 + delta;
        return {
          s: 3,
          pieces: [
            { cell: 1, u: a },
            { cell: 3, u: b },
          ],
        };
      };
      expect(verdict(d, "cube", pair(0)).ok).toBe(true);
      expect(verdict(d, "cube", pair(-2e-9)).ok).toBe(true);
      expect(verdict(d, "cube", pair(2e-9)).ok).toBe(false);
      expect(verdict(d, "cube", pair(5e-10), 1e-9).ok).toBe(true);
      expect(verdict(d, "cube", pair(5e-10), 0).ok).toBe(false);
    }
  });

  it("refuses malformed and out-of-scope input with reasons, and never throws on it", () => {
    const d = 6;
    const P = problem(d, "ball", 1);
    const ok = { s: 3, pieces: [{ cell: 0, u: zeros(5) }] };
    expect(validateNCube(P, ok).ok).toBe(true);
    const bad: unknown[] = [
      null,
      [],
      { s: NaN, pieces: ok.pieces },
      { s: -1, pieces: ok.pieces },
      { s: 3 },
      { s: 3, pieces: [] },
      { s: 3, pieces: [{ cell: 12, u: zeros(5) }] },
      { s: 3, pieces: [{ cell: 0, u: zeros(3) }] },
      { s: 3, pieces: [{ cell: 0, u: [0, 0, 0, 0, Infinity] }] },
      { s: 3, pieces: [{ cell: 0, u: zeros(5), q: [1, 0, 0, 0] }] },
      { s: 3, pieces: [{ cell: 0, u: [9, 0, 0, 0, 0] }] },
      { s: 3, pieces: [Object.assign([], { length: 1 })] },
    ];
    for (const c of bad) {
      const v = validateNCube(P, c);
      expect(v.ok).toBe(false);
      expect(v.reasons.length).toBeGreaterThan(0);
    }
    expect(() => validateNCube(P, ok, { tol: -1 })).toThrow(RangeError);
    expect(validateNCube({ polytope: "11cube" as never, piece: "ball", P: 1 }, ok).ok).toBe(false);
    // a sparse list: the hole is refused
    const holes = new Array(2);
    holes[0] = ok.pieces[0];
    expect(validateNCube(problem(d, "ball", 2), { s: 3, pieces: holes }).ok).toBe(false);
  });

  it("is capped: a check that needs more than its work refuses with a reason, quickly", () => {
    const d = 10;
    const n = d - 1;
    const a = 0.5 - 0.3536;
    const pieces = Array.from({ length: 12 }, (_, i) => ({ cell: i % 20, u: new Array(n).fill(a * (1 - i * 1e-3)) }));
    const v = validateNCube(problem(d, "ball", 12), { s: 1, pieces }, { maxWork: 40, maxReasons: 3 });
    expect(v.ok).toBe(false);
    expect(v.reasons.join(" ")).toMatch(/work limit|overlap/);
    const t0 = performance.now();
    const w = validateNCube(problem(d, "ball", 12), { s: 1, pieces });
    expect(performance.now() - t0).toBeLessThan(5000);
    expect(w.ok).toBe(false);
    expect(NCUBE_MAX_WORK).toBe(4_000_000);
  });
});

describe("the d-cube problems", () => {
  it("have ids of the boundary form, P up to 160, 384 and 400", () => {
    expect(parseBoundaryId("boundary/5cube/ball/160")).toEqual({ polytope: "5cube", piece: "ball", P: 160 });
    expect(parseBoundaryId("boundary/10cube/cube/400")).toEqual({ polytope: "10cube", piece: "cube", P: 400 });
    for (const bad of [
      "boundary/11cube/ball/1",
      "boundary/4cube/ball/1",
      "boundary/5cube/sphere/1",
      "boundary/5cube/ball/05",
    ])
      expect(typeof parseBoundaryId(bad)).toBe("string");
    expect(NCUBE_MAX_P).toEqual({
      "5cube": 160,
      "6cube": 384,
      "7cube": 400,
      "8cube": 400,
      "9cube": 400,
      "10cube": 400,
    });
    expect(boundaryDomain.bound({ fields: { polytope: "5cube", piece: "cube" }, P: 160 })!.value).toBe(2);
    expect(boundaryDomain.bound({ fields: { polytope: "5cube", piece: "cube" }, P: 161 })).toBeUndefined();
    expect(boundaryDomain.bound({ fields: { polytope: "9cube", piece: "ball" }, P: 400 })!.label).toMatch(
      /volume bound/,
    );
  });

  it("have the exact lower bounds: one ball 1/sqrt 2, one cube 1, a tiling's k", () => {
    for (let d = 5; d <= 10; d++) {
      expect(ncubeBound(d, "ball", 1)).toBe(boundaryBound("ball", 1));
      expect(ncubeBound(d, "ball", 1)).toBe(0.7071067811865475);
      expect(ncubeBound(d, "cube", 1)).toBe(1);
      expect(ncubeBound(d, "cube", 2 * d)).toBe(1);
      expect(ncubeBound(d, "cube", 2 * d + 1)).toBeGreaterThan(1);
      const n = d - 1;
      // P = 2d 2^n cubes need s >= 2 exactly (x^n <= 2^n in integers)
      if (2 * d * 2 ** n <= NCUBE_MAX_P[`${d}cube` as keyof typeof NCUBE_MAX_P])
        expect(ncubeBound(d, "cube", 2 * d * 2 ** n)).toBe(2);
      // a bound is below the packing that realises it
      const L = ncubeCellLayout(d, "ball", 4 * d, mulberry32(1));
      expect(ncubeBound(d, "ball", 4 * d)).toBeLessThanOrEqual(L.s);
    }
  });

  it("classify heavy checks by their estimated cost", () => {
    expect(boundaryIsHeavy({ polytope: "10cube", piece: "ball", P: 400 })).toBe(true);
    expect(boundaryIsHeavy({ polytope: "5cube", piece: "cube", P: 160 })).toBe(false);
    expect(boundaryIsHeavy({ polytope: "10cube", piece: "ball", P: 20 })).toBe(false);
  });

  it("answer the server's rule: P, the form, the bound, the packing", () => {
    const P = problem(5, "cube", 10);
    const L = packings.find((p) => p.id === "boundary/5cube/cube/10")!.config;
    expect(ncubeServerRule(P, L, 1).valid).toBe(true);
    expect(ncubeServerRule({ ...P, P: 161 }, L, 1).reasons[0]).toMatch(/above the largest P/);
    expect(ncubeServerRule(P, { ...L, s: 0.99 }, 1).reasons[0]).toMatch(/below the lower bound/);
    expect(
      ncubeServerRule(P, { s: 1.01, pieces: L.pieces.map((q) => ({ ...q, q: [1, 0, 0, 0] })) }, 1).reasons[0],
    ).toMatch(/no rotation q/);
    expect(
      ncubeServerRule(P, { s: 1.01, pieces: L.pieces.map((q) => ({ ...q, u: [0, 0, 0] })) }, 1).reasons[0],
    ).toMatch(/4 numbers/);
  });
});

describe("the d-cube kernel", () => {
  const kernel = (d: number, piece: Piece, P: number) => boundaryKernel(problem(d, piece, P));

  it("passes the contract (small problems; the larger take minutes to check)", () => {
    for (const [d, piece, P] of [
      [5, "ball", 1],
      [5, "ball", 6],
      [5, "cube", 1],
      [5, "cube", 12],
      [6, "ball", 3],
      [6, "cube", 7],
      [7, "ball", 2],
      [8, "cube", 5],
      [10, "ball", 1],
      [10, "ball", 2],
      [10, "cube", 9],
    ] as const)
      expect(
        checkKernel(() => kernel(d, piece, P), { seed: P, states: 2 }),
        `${d}cube ${piece} ${P}`,
      ).toEqual([]);
  }, 300_000);

  it("seeds valid structured layouts (cell layouts; for cubes and valid balls the ring cut), each with no energy", () => {
    const rng = mulberry32(4);
    for (const d of [5, 6, 7, 10])
      for (const piece of ["ball", "cube"] as const)
        for (const P of [3, 2 * d, 3 * d + 1, 40]) {
          const k = kernel(d, piece, P);
          const seeds = k.seeds!(rng, 4);
          expect(seeds).toHaveLength(4);
          for (const q of seeds) {
            expect(k.energy(q), `${d} ${piece} ${P}`).toBeLessThan(1e-20);
            expect(k.validate(k.toConfig(q)).reasons).toEqual([]);
          }
        }
    // one ball per cell at s = 1; cubes: one per cell at 1 + 1e-6
    expect(kernel(6, "ball", 12).seeds!(rng, 2).map((q) => q.size)).toContain(1);
    expect(kernel(6, "cube", 12).seeds!(rng, 2).map((q) => q.size)).toContain(1 + 1e-6);
    // the ring cut straddles ridges: cubes at s = m / 4 round the ring (P = 20 -> s = 5/4 at a cross grid of 1... or a cell layout if smaller)
    const ring = ncubeRingLayout(5, "cube", 12, mulberry32(3))!;
    expect(ring.s).toBeLessThanOrEqual(3 + 1e-6);
    const cfg = {
      s: ring.s,
      pieces: Array.from({ length: 12 }, (_, i) => ({
        cell: ring.st.cells[i],
        u: Array.from(ring.st.U.slice(4 * i, 4 * i + 4)),
      })),
    };
    expect(verdict(5, "cube", cfg).reasons).toEqual([]);
  });

  it("finishes a squeezed start into a valid configuration and returns the same from the same seed", () => {
    for (const [d, piece, P] of [
      [5, "ball", 4],
      [5, "cube", 6],
      [7, "ball", 3],
    ] as const) {
      const k = kernel(d, piece, P);
      const s = k.randomStart(mulberry32(100 + P));
      squeeze(k, s);
      const f = k.finish(s)!;
      expect(f, `${d} ${piece} ${P}`).not.toBeNull();
      expect(f.score).toBeGreaterThanOrEqual(k.bound! * (1 - 1e-12));
      expect(verdict(d, piece, f.config as Cfg).reasons).toEqual([]);
      expect(k.validate(f.config).valid).toBe(true);
    }
  });

  it("has the gradient of its weighted energy, in every displacement and the size, across ridges", () => {
    for (const [d, piece, P] of [
      [5, "ball", 4],
      [5, "cube", 4],
      [6, "ball", 3],
    ] as const) {
      const k = kernel(d, piece, P);
      const rng = mulberry32(P + d);
      const n = (d - 1) * P;
      for (let t = 0; t < 4; t++) {
        const s = k.randomStart(rng);
        k.setSize(s, s.size * (0.8 + 0.2 * rng()));
        for (let i = 0; i < P; i++) k.jiggle(s, i, (1.5 * s.size) / k.steps[0], rng);
        const w = Float64Array.from({ length: P }, () => 0.5 + rng());
        const g = new Float64Array(n + 1);
        if (!(k.energy(s, g, w) > 1e-8)) continue;
        for (let i = 0; i <= n; i++) {
          const h = 1e-7 * (i === n ? s.size : 1);
          const p = k.clone(s);
          const m = k.clone(s);
          if (i === n) {
            p.size += h;
            m.size -= h;
          } else {
            p.x[i] += h;
            m.x[i] -= h;
          }
          const fd = (k.energy(p, undefined, w) - k.energy(m, undefined, w)) / (2 * h);
          expect(Math.abs(fd - g[i]), `${d} ${piece} i=${i}`).toBeLessThan(1e-4 * Math.max(1, Math.abs(fd)));
        }
      }
    }
  });

  it("refuses a P above the problem set's and an unknown container", () => {
    expect(() => kernel(5, "ball", 161)).toThrow(/above the largest/);
    expect(() => boundaryKernel({ polytope: "11cube" as never, piece: "ball", P: 1 })).toThrow();
  });
});

describe("the loop lemma (BOUNDARY.md 4B.3)", () => {
  it("never finds a loop shorter than 1 once s is at least 1/sqrt 2, so the validator skips the loop search there", () => {
    for (const d of [5, 6, 7]) {
      const n = d - 1;
      const rnd = mulberry32(40 + d);
      let tried = 0;
      for (let t = 0; t < 120; t++) {
        const s = Math.SQRT1_2 + 1.2 * rnd() ** 2;
        const h = s / 2;
        const tight = rnd() < 0.5;
        const u = Array.from({ length: n }, () =>
          tight ? (rnd() < 0.5 ? -1 : 1) * Math.max(0, h - 0.3536) * (0.9 + 0.1 * rnd()) : (rnd() * 2 - 1) * h,
        );
        if (ncubeClearance(Float64Array.from(u), 0, n, h) < 0.5) continue;
        tried++;
        const c = Math.floor(rnd() * 2 * d);
        let best = Infinity;
        for (const r of ncubeRealisedWords(
          d,
          s,
          c,
          Float64Array.from(u),
          0,
          c,
          Float64Array.from(u),
          0,
          ncubeBudget(2e6),
          {
            excludeEmpty: true,
          },
        ))
          if (r.len > 1e-12) best = Math.min(best, r.len);
        expect(best, `d=${d} s=${s}`).toBeGreaterThanOrEqual(1);
      }
      expect(tried).toBeGreaterThan(15);
    }
  });

  it("still rejects a self-overlapping ball below 1/sqrt 2 (the loop search runs there)", () => {
    // the loop check needs s sqrt 2 < 1 to be able to fire; at such an s the clearance rule (a ball of diameter 1 clear of the faces) already fails
    // for every u, so the search is only a safeguard: the verdict is invalid either way
    const v = validateBoundary(problem(5, "ball", 1), { s: 0.7, pieces: [{ cell: 0, u: [0, 0, 0, 0] }] });
    expect(v.ok).toBe(false);
  });
});
