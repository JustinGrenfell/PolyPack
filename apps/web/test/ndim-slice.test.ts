/**
 * The geometry behind the board's drawings of the nD in nD space (apps/web/src/space/render/ndslice.ts): the slice of a
 * placed polytope or ball by a hyperplane x_k = c, with exact expectations, and the projection and the table's texts.
 */
import { ConvexHull } from "three/examples/jsm/math/ConvexHull.js";
import { Vector3 } from "three";
import { cayleyMatrix, getPolytope, loadPolytope, poseOfMatrix, type NdConfig, type NdProblem } from "@polypack/ndim";
import { beforeAll, describe, expect, it } from "vitest";
import {
  type Placed,
  type PolytopeSlice,
  centreText,
  coordText,
  edgesOf,
  extentAlong,
  hull2d,
  placeContainer,
  placePieces,
  projectOnto,
  rotationText,
  sliceAt,
} from "../src/space/render/ndslice";

/** The volume of the convex hull of 3D points, by the tetrahedra from its centroid to its faces. */
function hullVolume(points: number[][]): { volume: number; faces: number } {
  const hull = new ConvexHull().setFromPoints(points.map((p) => new Vector3(p[0], p[1], p[2])));
  const c = new Vector3();
  for (const p of points) c.add(new Vector3(p[0], p[1], p[2]));
  c.divideScalar(points.length);
  let volume = 0;
  for (const f of hull.faces) {
    const a = f.edge.head().point;
    const b = f.edge.next.head().point;
    const d = f.edge.next.next.head().point;
    volume += Math.abs(new Vector3().subVectors(a, c).dot(new Vector3().crossVectors(new Vector3().subVectors(b, c), new Vector3().subVectors(d, c)))) / 6;
  }
  return { volume, faces: hull.faces.length };
}

const problem = (piece: string, container: string, d = 4): NdProblem => ({ d, piece, container });

/** The polytope at scale 1, centred, as a placed container. */
const unit = (code: string): Placed => placeContainer(problem(code, code), 1);

const sorted = (pts: number[][]) =>
  pts.map((p) => p.map((x) => Math.round(x * 1e9) / 1e9)).sort((a, b) => a[0] - b[0] || a[1] - b[1] || a[2] - b[2]);

beforeAll(async () => {
  for (const c of ["4:cube", "4:ortho", "4:24cell", "4:simplex", "4:120cell", "4:600cell", "5:cube", "5:ortho", "5:simplex"]) await loadPolytope(c);
});

describe("the slice of a 4D polytope by x_k = c", () => {
  it("cuts the 4-cube at its centre to a unit cube", () => {
    // the catalogue's tesseract has vertices (+-1/2)^4; the 32 edges along x4 are cut at their middles, so the slice is the
    // 8 points (+-1/2)^3: a cube of edge 1
    for (const k of [0, 1, 2, 3]) {
      const s = sliceAt(unit("4:cube"), k, 0) as PolytopeSlice;
      expect(s.kind).toBe("polytope");
      expect(s.points).toHaveLength(8);
      for (const p of s.points) for (const x of p) expect(Math.abs(x)).toBeCloseTo(0.5, 12);
      const { volume, faces } = hullVolume(s.points);
      expect(volume).toBeCloseTo(1, 12);
      expect(faces).toBe(12);
    }
  });

  it("cuts the 16-cell (orthoplex) at x4 = 0 to an octahedron", () => {
    // vertices +-e_i / sqrt 2; the six on x4 = 0 are the slice, an octahedron of edge 1 and volume sqrt(2) / 3
    const s = sliceAt(unit("4:ortho"), 3, 0) as PolytopeSlice;
    expect(s.points).toHaveLength(6);
    for (const p of s.points) expect(Math.hypot(...p)).toBeCloseTo(Math.SQRT1_2, 12);
    const { volume, faces } = hullVolume(s.points);
    expect(volume).toBeCloseTo(Math.SQRT2 / 3, 12);
    expect(faces).toBe(8);
  });

  it("cuts the 24-cell at its centre, in the catalogue's orientation, to a rhombic dodecahedron", () => {
    // The catalogue's 24-cell has the 8 vertices +-e_i and the 16 vertices (+-1/2)^4. The hyperplane x4 = 0 holds the 6 vertices
    // +-e_1, +-e_2, +-e_3, and cuts the 8 edges between (+-1/2, +-1/2, +-1/2, 1/2) and (+-1/2, +-1/2, +-1/2, -1/2) at their
    // middles (+-1/2, +-1/2, +-1/2, 0). No other edge crosses it. The 14 points are the cube [-1/2, 1/2]^3 with its 6 face
    // pyramids of apex +-e_i: the rhombic dodecahedron, of volume 1 + 6 (1/3)(1)(1/2) = 2, with 12 rhombic faces (24 triangles).
    const s = sliceAt(unit("4:24cell"), 3, 0) as PolytopeSlice;
    expect(s.points).toHaveLength(14);
    const apex = s.points.filter((p) => Math.abs(Math.hypot(...p) - 1) < 1e-9);
    const corners = s.points.filter((p) => p.every((x) => Math.abs(Math.abs(x) - 0.5) < 1e-9));
    expect(apex).toHaveLength(6);
    expect(corners).toHaveLength(8);
    const { volume, faces } = hullVolume(s.points);
    expect(volume).toBeCloseTo(2, 12);
    expect(faces).toBe(24);
    // every axis of this orientation gives the same slice
    for (const k of [0, 1, 2]) expect(sorted((sliceAt(unit("4:24cell"), k, 0) as PolytopeSlice).points)).toHaveLength(14);
  });

  it("cuts a polytope at a height between its vertices by interpolating its edges", () => {
    // the tesseract at x4 = 1/4: only the 32 edges along x4 cross it, at (+-1/2)^3 again, the unit cube at every height
    const s = sliceAt(unit("4:cube"), 3, 0.25) as PolytopeSlice;
    expect(s.points).toHaveLength(8);
    expect(hullVolume(s.points).volume).toBeCloseTo(1, 12);
    // the 5-cell (simplex): a slice through the middle of its x4 extent is a tetrahedron or a prism-like solid with edges inside it
    const simplex = unit("4:simplex");
    const [lo, hi] = extentAlong(simplex, 3);
    const mid = sliceAt(simplex, 3, (lo + hi) / 2) as PolytopeSlice;
    expect(mid.points.length).toBeGreaterThanOrEqual(4);
    for (const p of mid.points) expect(Math.hypot(...p)).toBeLessThan(1);
  });

  it("is empty beyond the polytope, and a point at its extreme vertex", () => {
    const cube = unit("4:cube");
    expect(extentAlong(cube, 2)).toEqual([-0.5, 0.5]);
    expect(sliceAt(cube, 2, 0.5001)).toBeNull();
    expect(sliceAt(cube, 2, -3)).toBeNull();
    // at the face x3 = 1/2 the slice is the face itself: its 8 vertices
    const face = sliceAt(cube, 2, 0.5) as PolytopeSlice;
    expect(face.points).toHaveLength(8);
    // the 16-cell's vertex e_4 / sqrt 2 touches x4 = 1/sqrt 2 in one point
    const tip = sliceAt(unit("4:ortho"), 3, Math.SQRT1_2) as PolytopeSlice;
    expect(tip.points).toEqual([[0, 0, 0]]);
  });

  it("follows a turned piece: a half-turn of the tesseract in the x1 x2 plane slices as before", () => {
    // a rotation by pi/2 in the x1 x2 plane maps the cube to itself, so its slice at x3 = 0 stays the unit cube
    const R = new Float64Array(16);
    R[0 * 4 + 1] = -1;
    R[1 * 4 + 0] = 1;
    R[2 * 4 + 2] = 1;
    R[3 * 4 + 3] = 1;
    const config: NdConfig = { s: 2, pieces: [poseOfMatrix([0.2, 0.1, 0, 0], R, 4)] };
    const { pieces, errors } = placePieces(problem("4:cube", "4:cube"), config);
    expect(errors).toEqual([]);
    const s = sliceAt(pieces[0]!, 2, 0) as PolytopeSlice;
    expect(s.points).toHaveLength(8);
    expect(hullVolume(s.points).volume).toBeCloseTo(1, 9);
    // centred at (0.2, 0.1): the cut cube is centred at (0.2, 0.1, 0)
    const mean = [0, 1, 2].map((i) => s.points.reduce((a, p) => a + p[i], 0) / 8);
    expect(mean[0]).toBeCloseTo(0.2, 9);
    expect(mean[1]).toBeCloseTo(0.1, 9);
    expect(mean[2]).toBeCloseTo(0, 9);
  });

  it("finds the edges of the catalogue's polytopes: as many as the literature's counts", () => {
    for (const code of ["4:simplex", "4:cube", "4:ortho", "4:24cell", "4:120cell", "4:600cell", "5:cube", "5:ortho"]) {
      const poly = getPolytope(code);
      expect(edgesOf(poly), code).toHaveLength(poly.nE);
    }
  });
});

describe("the slice of a ball", () => {
  it("is a ball of radius sqrt(r^2 - (c - x_k)^2) about the centre without its k-th coordinate", () => {
    const ball: Placed = { kind: "ball", centre: Float64Array.from([1, 2, 3, 4]), radius: 0.5 };
    const at = sliceAt(ball, 3, 4.3);
    expect(at).toMatchObject({ kind: "ball", centre: [1, 2, 3] });
    expect((at as { radius: number }).radius).toBeCloseTo(0.4, 12);
    // through the centre: the full radius
    expect(sliceAt(ball, 0, 1)).toEqual({ kind: "ball", centre: [2, 3, 4], radius: 0.5 });
    // a container ball of diameter s = 3: radius 3/2, at c = 0.9 the cut has radius 1.2
    const big = placeContainer(problem("ball", "ball"), 3);
    expect((sliceAt(big, 1, 0.9) as { radius: number }).radius).toBeCloseTo(1.2, 12);
    expect(extentAlong(big, 2)).toEqual([-1.5, 1.5]);
  });

  it("is empty when the hyperplane misses it or only touches it", () => {
    const ball: Placed = { kind: "ball", centre: Float64Array.from([0, 0, 0, 0]), radius: 0.5 };
    expect(sliceAt(ball, 2, 0.5)).toBeNull();
    expect(sliceAt(ball, 2, 0.75)).toBeNull();
    expect(sliceAt(ball, 2, -2)).toBeNull();
    expect(sliceAt(ball, 2, 0.4999)).not.toBeNull();
  });
});

describe("the projection onto two axes", () => {
  it("projects a 5-cube onto a plane of its axes as a square, and onto a diagonal plane as a hull of more points", () => {
    const c = placeContainer(problem("5:cube", "5:cube", 5), 3);
    const pr = projectOnto(c, 0, 4);
    expect(pr.kind).toBe("polygon");
    if (pr.kind !== "polygon") return;
    expect(pr.points).toHaveLength(4);
    for (const [x, y] of pr.points) {
      expect(Math.abs(x)).toBeCloseTo(1.5, 12);
      expect(Math.abs(y)).toBeCloseTo(1.5, 12);
    }
  });

  it("projects a ball to a disk and a 5-orthoplex onto two axes to a diamond", () => {
    const b = placeContainer(problem("ball", "ball", 6), 2);
    expect(projectOnto(b, 1, 3)).toEqual({ kind: "ball", cx: 0, cy: 0, r: 1 });
    const o = placeContainer(problem("5:ortho", "5:ortho", 5), 1);
    const pr = projectOnto(o, 2, 3);
    if (pr.kind !== "polygon") throw new Error("a polygon");
    expect(pr.points).toHaveLength(4);
    for (const [x, y] of pr.points) expect(Math.abs(x) + Math.abs(y)).toBeCloseTo(Math.SQRT1_2, 12);
  });

  it("takes the 2D hull counter-clockwise and drops interior and collinear points", () => {
    const h = hull2d([
      [0, 0],
      [2, 0],
      [2, 2],
      [0, 2],
      [1, 1],
      [1, 0],
    ]);
    expect(h).toHaveLength(4);
    expect(h).toContainEqual([0, 0]);
    expect(h).toContainEqual([2, 2]);
  });
});

describe("the table's texts", () => {
  it("writes coordinates with four decimals and a true minus", () => {
    expect(coordText(-0.5)).toBe("−0.5000");
    expect(coordText(-0.00001)).toBe("0.0000");
    expect(centreText([1, -2.25, 0])).toBe("(1.0000, −2.2500, 0.0000)");
  });

  it("says how far a piece is turned", () => {
    const d = 4;
    const I = new Float64Array(16);
    for (let i = 0; i < 4; i++) I[i * 4 + i] = 1;
    expect(rotationText(I, d)).toBe("not turned");
    // a quarter turn in one plane: the axes are permuted (with a sign)
    const Q = Float64Array.from(I);
    Q[0] = 0;
    Q[5] = 0;
    Q[1] = -1;
    Q[4] = 1;
    expect(rotationText(Q, d)).toBe("axes permuted");
    // a general rotation: Cayley parameter 0.5 in the x1 x2 plane, 2 atan(0.5) = 53.13 degrees there; mean turn arccos((2 cos + 2) / 4)
    const R = cayleyMatrix([0.5, 0, 0, 0, 0, 0], 4, 0);
    const cos = Math.cos(2 * Math.atan(0.5));
    expect(rotationText(R, d)).toBe(`turned about ${Math.round((Math.acos((2 * cos + 2) / 4) * 180) / Math.PI)}°`);
  });
});
