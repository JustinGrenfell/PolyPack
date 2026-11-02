/**
 * The board's drawing of a 4D in 4D to 10D in 10D packing (ndim/<d>/<piece>/<container>/<P>, docs/spaces/NDIM.md).
 *
 *  - 4D in 4D: a 3D slice. The container and every piece are cut by the hyperplane x_k = c; the user picks the axis k and
 *    moves c over the container's extent. A cut polytope is a convex polyhedron (ndslice.ts), drawn with three.js in the
 *    look of the polyhedra drawing; a cut ball is a ball.
 *  - 5D to 10D: a chart (an orthographic projection of the container and the pieces onto two axes the user picks) and a
 *    table of the pieces.
 *
 * Both load their polytope tables with loadPolytope (the lazy ones come as chunks of the ndim package) and three.js only
 * when the first 4D drawing is made. The geometry is light enough for the page's thread; validation and search run in the
 * worker (search.ts).
 */
import type * as T from "three";
import { BALL, type NdConfig, type NdProblem } from "@polypack/ndim";
import { ndCode, ndPairTitle, ndShapeTitle } from "@polypack/ndspace/catalogue";
import type { SpaceConfig } from "@polypack/spaces";
import type { SearchProblem } from "@polypack/spaces/search";
import type { SpaceView } from "../view";
import { cssVar, isDark, onDrag, pieceFill } from "./common";
import {
  type Placed,
  centreText,
  extentAlong,
  loadShapes,
  placeContainer,
  placePieces,
  projectOnto,
  rotationText,
  sliceAt,
} from "./ndslice";
import { chartColours, dataTable, h, labelledSelect, svg } from "./nd-ui";

type NdimProblem = Extract<SearchProblem, { space: "ndim" }>;
type Three = typeof import("three");

/** The 4D drawing's starting turn: 25 degrees about the vertical axis, tilted 20 degrees (as polyhedra.ts). */
const START = { yaw: (25 * Math.PI) / 180, pitch: (-20 * Math.PI) / 180 };

const isConfig = (c: unknown): c is NdConfig =>
  typeof c === "object" && c !== null && typeof (c as NdConfig).s === "number" && Array.isArray((c as NdConfig).pieces);

/** "12 hypercubes in an orthoplex of scale 2.5". */
const describeConfig = (p: NdimProblem, s: number): string =>
  `${p.P} ${ndPairTitle(p.piece, p.container).toLowerCase()} of scale ${s.toFixed(6)}`;

export function mountNdimView(host: HTMLElement, p: NdimProblem, label = "Drawing"): SpaceView {
  const d = p.d;
  // the catalogue's codes ("4:cube") of the problem's shapes ("cube"), which the geometry reads
  const np: NdProblem = { d, piece: ndCode(d, p.piece), container: ndCode(d, p.container) };
  const root = h("div", { class: `space-view space-view-nd space-view-nd${d === 4 ? "-slice" : "-chart"}` });
  const tools = h("div", { class: "toolbar space-tools nd-tools", role: "toolbar", "aria-label": label });
  const status = h("p", { class: "note nd-status", role: "status" }, "Loading the polytopes…");
  const box = h("div", { class: "space-draw nd-draw" });
  const tableHost = h("div", { class: "nd-pieces" });
  root.append(tools, status, box, tableHost);
  host.replaceChildren(root);

  let config: NdConfig | null = null;
  let container: Placed | null = null;
  let pieces: (Placed | null)[] = [];
  let ready = false;
  let gone = false;
  let wanted: SpaceConfig | null | undefined;

  const axes = Array.from({ length: d }, (_, k): [string, string] => [String(k), `x${k + 1}`]);
  // the slice axis and position (4D), or the two axes of the chart
  let axis = d - 1;
  let at = 0;
  let ax = 0;
  let ay = 1;

  const slice = d === 4 ? sliceDrawing() : null;
  const chart = d === 4 ? null : chartDrawing();

  // ---------- controls ----------
  let readout: HTMLElement | null = null;
  let slider: HTMLInputElement | null = null;
  let hint: HTMLElement | null = null;
  if (d === 4) {
    const axisField = labelledSelect(
      "Slice at",
      axes,
      String(axis),
      (v) => {
        axis = Number(v);
        at = 0;
        syncSlider();
        draw();
      },
      { "data-nd": "axis" },
    );
    slider = h("input", {
      type: "range",
      class: "nd-slider",
      "aria-label": "Position of the slice along the axis",
      "data-nd": "position",
      min: "-1",
      max: "1",
      step: "any",
      value: "0",
    });
    slider.addEventListener("input", () => {
      at = Number(slider!.value);
      readoutText();
      draw();
    });
    readout = h("output", { class: "nd-readout", "data-nd": "readout" });
    const centre = h("button", { type: "button" }, "Centre");
    centre.addEventListener("click", () => {
      at = 0;
      syncSlider();
      draw();
    });
    const reset = h("button", { type: "button" }, "Reset view");
    reset.addEventListener("click", () => slice!.reset());
    hint = h("span", { class: "space-hint" }, "Drag to turn");
    tools.append(axisField.box, slider, readout, centre, reset, hint);
  } else {
    const fx = labelledSelect("Horizontal axis", axes, String(ax), (v) => {
      ax = Number(v);
      if (ax === ay) {
        ay = (ax + 1) % d;
        fy.select.value = String(ay);
      }
      draw();
    }, { "data-nd": "axis-x" });
    const fy = labelledSelect("Vertical axis", axes, String(ay), (v) => {
      ay = Number(v);
      if (ax === ay) {
        ax = (ay + 1) % d;
        fx.select.value = String(ax);
      }
      draw();
    }, { "data-nd": "axis-y" });
    tools.append(fx.box, fy.box);
  }

  function readoutText(): void {
    if (readout) readout.textContent = `x${axis + 1} = ${at.toFixed(3).replace(/^-/, "−")}`;
  }

  /** The slider covers the container's extent along the axis; c stays inside it. */
  function syncSlider(): void {
    if (!slider) return;
    if (container) {
      const [lo, hi] = extentAlong(container, axis);
      slider.min = String(lo);
      slider.max = String(hi);
      slider.step = String((hi - lo) / 400);
      at = Math.min(hi, Math.max(lo, at));
    }
    slider.value = String(at);
    readoutText();
  }

  // ---------- the pieces' table ----------
  function renderTable(): void {
    if (!config) {
      tableHost.replaceChildren();
      return;
    }
    const rows = pieces.map((q, i) => [
      String(i + 1),
      ndShapeTitle(p.piece),
      q ? centreText(q.centre) : "unreadable",
      q?.kind === "polytope" ? rotationText(q.R, d) : p.piece === BALL ? "none (a ball)" : "—",
    ]);
    const table = dataTable(
      `Pieces of the packing (${p.P}), centre coordinates x1 to x${d}`,
      ["#", "Piece", "Centre", "Turn"],
      rows,
      [0],
    );
    if (d === 4) {
      const details = h("details", { class: "nd-details" }, h("summary", {}, "Pieces"), table);
      tableHost.replaceChildren(details);
    } else tableHost.replaceChildren(table);
  }

  // ---------- drawing ----------
  function draw(): void {
    if (gone || !ready) return;
    if (!config || !container) {
      status.textContent = "No packing to draw.";
      status.hidden = false;
      box.style.visibility = "hidden";
      return;
    }
    status.hidden = true;
    box.style.visibility = "";
    box.setAttribute("aria-label", describe());
    if (slice) slice.draw(container, pieces, axis, at, config.s);
    else chart!.draw(container, pieces, ax, ay);
  }

  const describe = (): string =>
    config
      ? d === 4
        ? `${describeConfig(p, config.s)}, cut by x${axis + 1} = ${at.toFixed(3)}, drawn in 3D.`
        : `${describeConfig(p, config.s)}, projected onto x${ax + 1} and x${ay + 1}.`
      : "No packing to draw.";

  function place(): void {
    if (!isConfig(wanted)) {
      config = null;
      container = null;
      pieces = [];
    } else {
      config = wanted;
      container = placeContainer(np, config.s);
      const placed = placePieces(np, config);
      pieces = placed.pieces;
    }
    syncSlider();
    renderTable();
    draw();
  }

  box.setAttribute("role", "img");
  if (slice) {
    box.tabIndex = 0;
    box.setAttribute("role", "application");
    box.append(slice.element);
    slice.attach(box);
    box.addEventListener("keydown", (e) => {
      const k = e.shiftKey ? 0.15 : 0.05;
      const move: Record<string, [number, number]> = {
        ArrowLeft: [-k, 0],
        ArrowRight: [k, 0],
        ArrowUp: [0, -k],
        ArrowDown: [0, k],
      };
      if (move[e.key]) {
        e.preventDefault();
        slice.turn(...move[e.key]);
      } else if (e.key === "Home" || e.key === "0") {
        e.preventDefault();
        slice.reset();
      }
    });
    hint!.hidden = true;
  } else box.append(chart!.element);

  const resize = new ResizeObserver(() => draw());
  resize.observe(box);
  const scheme = matchMedia("(prefers-color-scheme: dark)");
  scheme.addEventListener("change", draw);

  loadShapes(np).then(
    () => {
      if (gone) return;
      ready = true;
      if (wanted !== undefined) place();
      else status.textContent = "No packing to draw.";
    },
    (error: unknown) => {
      status.textContent = `The polytopes could not be loaded: ${error instanceof Error ? error.message : String(error)}.`;
    },
  );

  return {
    show(c) {
      wanted = c;
      if (ready) place();
      if (hint) hint.hidden = c === null;
    },
    // the same shapes in another P keep the axis, the slice and the turn the user chose
    canShow: (q) => q.space === "ndim" && q.d === p.d && q.piece === p.piece && q.container === p.container,
    destroy() {
      gone = true;
      resize.disconnect();
      scheme.removeEventListener("change", draw);
      slice?.destroy();
      root.remove();
    },
  };
}

// ---------- the chart of 5D to 10D ----------

function chartDrawing() {
  const element = svg("svg", { class: "nd-chart", role: "presentation", preserveAspectRatio: "xMidYMid meet" });
  return {
    element,
    draw(container: Placed, pieces: (Placed | null)[], i: number, j: number) {
      const col = chartColours();
      const frame = projectOnto(container, i, j);
      // the plot's box: the container's projection, with a margin
      let x0: number, x1: number, y0: number, y1: number;
      if (frame.kind === "ball") {
        [x0, x1, y0, y1] = [frame.cx - frame.r, frame.cx + frame.r, frame.cy - frame.r, frame.cy + frame.r];
      } else {
        const xs = frame.points.map((q) => q[0]);
        const ys = frame.points.map((q) => q[1]);
        [x0, x1, y0, y1] = [Math.min(...xs), Math.max(...xs), Math.min(...ys), Math.max(...ys)];
      }
      const span = Math.max(x1 - x0, y1 - y0) || 1;
      const m = 0.06 * span;
      // y points up on the chart, down in SVG
      element.setAttribute("viewBox", `${x0 - m} ${-(y1 + m)} ${x1 - x0 + 2 * m} ${y1 - y0 + 2 * m}`);
      const stroke = { "vector-effect": "non-scaling-stroke", "stroke-linejoin": "round" } as const;
      const shape = (pr: ReturnType<typeof projectOnto>, attrs: Record<string, string | number>) =>
        pr.kind === "ball"
          ? svg("circle", { cx: pr.cx, cy: -pr.cy, r: pr.r, ...stroke, ...attrs })
          : svg("polygon", { points: pr.points.map((q) => `${q[0]},${-q[1]}`).join(" "), ...stroke, ...attrs });
      const parts: SVGElement[] = [shape(frame, { fill: col.face, "fill-opacity": 0.5, stroke: col.edge, "stroke-width": 1.5 })];
      const fontSize = 0.035 * span;
      pieces.forEach((q, n) => {
        if (!q) return;
        const pr = projectOnto(q, i, j);
        const g = svg("g", { "data-piece": n + 1 });
        g.append(shape(pr, { fill: pieceFill(n, 0.5), stroke: col.piece, "stroke-width": 1 }));
        g.append(
          svg(
            "text",
            {
              x: q.centre[i],
              y: -q.centre[j],
              "font-size": fontSize,
              "text-anchor": "middle",
              "dominant-baseline": "central",
              fill: col.text,
            },
            String(n + 1),
          ),
        );
        parts.push(g);
      });
      element.replaceChildren(...parts);
    },
  };
}

// ---------- the slice of 4D in 4D ----------

function sliceDrawing() {
  const cv = document.createElement("canvas");
  let THREE: Three | null = null;
  let Hull: typeof import("three/examples/jsm/geometries/ConvexGeometry.js").ConvexGeometry | null = null;
  let failed = false;
  let destroyed = false;
  let lost = false;
  let gl: T.WebGLRenderer | null = null;
  let scene: T.Scene | null = null;
  let camera: T.PerspectiveCamera | null = null;
  let root: T.Group | null = null;
  let built: T.Object3D[] = [];
  let yaw = START.yaw;
  let pitch = START.pitch;
  let last: { container: Placed; pieces: (Placed | null)[]; axis: number; at: number; s: number } | null = null;
  let off: (() => void) | null = null;

  const dispose = () => {
    for (const o of built) {
      root?.remove(o);
      (o as T.Mesh).geometry?.dispose();
    }
    built = [];
  };

  function setup(T3: Three) {
    try {
      gl = new T3.WebGLRenderer({ canvas: cv, antialias: true, alpha: true });
    } catch {
      failed = true;
      return;
    }
    gl.setClearColor(0x000000, 0);
    scene = new T3.Scene();
    camera = new T3.PerspectiveCamera(30, 1, 0.01, 1e4);
    camera.add(new T3.AmbientLight(0xffffff, 1.4));
    const sun = new T3.DirectionalLight(0xffffff, 2.2);
    sun.position.set(-0.45, 0.6, 0.66);
    camera.add(sun);
    scene.add(camera);
    root = new T3.Group();
    scene.add(root);
  }

  /** The convex hull of 3D points as a geometry with flat faces, or null when they do not span a solid. */
  function hullGeometry(points: number[][]): T.BufferGeometry | null {
    if (!THREE || !Hull || points.length < 4) return null;
    try {
      const g = new Hull(points.map((q) => new THREE!.Vector3(q[0], q[1], q[2])));
      return g.getAttribute("position").count >= 12 ? g : (g.dispose(), null);
    } catch {
      // coplanar or collinear points: a touch, not a solid
      return null;
    }
  }

  function build() {
    if (!THREE || !root || !last) return;
    const T3 = THREE;
    dispose();
    const dark = isDark();
    const edgeColour = new T3.Color(cssVar("--space-edge") || (dark ? "#c9cbd0" : "#3b3b3f"));
    const faceColour = new T3.Color(cssVar("--space-face") || (dark ? "#2a2c31" : "#ece8df"));
    const edgeMaterial = new T3.LineBasicMaterial({ color: edgeColour, transparent: true, opacity: dark ? 0.55 : 0.5 });
    const boxEdge = new T3.LineBasicMaterial({ color: edgeColour, transparent: true, opacity: 0.9 });
    const boxFace = new T3.MeshBasicMaterial({
      color: faceColour,
      transparent: true,
      opacity: dark ? 0.1 : 0.16,
      depthWrite: false,
      side: T3.DoubleSide,
    });
    const add = (o: T.Object3D) => {
      root!.add(o);
      built.push(o);
    };
    const { container, pieces, axis, at } = last;
    // the container's cut: faint faces and its edges
    const cut = sliceAt(container, axis, at);
    if (cut?.kind === "ball") {
      const g = new T3.SphereGeometry(cut.radius, 48, 24);
      const m = new T3.Mesh(g, boxFace);
      m.position.set(cut.centre[0], cut.centre[1], cut.centre[2]);
      m.renderOrder = 2;
      add(m);
      for (let a = 0; a < 3; a++) {
        const ring = new T3.BufferGeometry().setFromPoints(
          Array.from({ length: 97 }, (_, k) => {
            const t = (2 * Math.PI * k) / 96;
            const u = [cut.radius * Math.cos(t), cut.radius * Math.sin(t), 0];
            const v = a === 0 ? [u[0], u[1], 0] : a === 1 ? [u[0], 0, u[1]] : [0, u[0], u[1]];
            return new T3.Vector3(cut.centre[0] + v[0], cut.centre[1] + v[1], cut.centre[2] + v[2]);
          }),
        );
        add(new T3.Line(ring, boxEdge));
      }
    } else if (cut?.kind === "polytope") {
      const g = hullGeometry(cut.points);
      if (g) {
        const m = new T3.Mesh(g, boxFace);
        m.renderOrder = 2;
        add(m);
        add(new T3.LineSegments(new T3.EdgesGeometry(g, 1), boxEdge));
      }
    }
    pieces.forEach((q, i) => {
      if (!q) return;
      const s = sliceAt(q, axis, at);
      if (!s) return;
      const colour = new T3.Color().setHSL(((i * 137.508) % 360) / 360, dark ? 0.45 : 0.55, dark ? 0.52 : 0.62, T3.SRGBColorSpace);
      const material = new T3.MeshLambertMaterial({
        color: colour,
        polygonOffset: true,
        polygonOffsetFactor: 1,
        polygonOffsetUnits: 1,
      });
      if (s.kind === "ball") {
        const m = new T3.Mesh(new T3.SphereGeometry(s.radius, 32, 16), material);
        m.position.set(s.centre[0], s.centre[1], s.centre[2]);
        add(m);
        return;
      }
      const g = hullGeometry(s.points);
      if (g) {
        add(new T3.Mesh(g, material));
        add(new T3.LineSegments(new T3.EdgesGeometry(g, 1), edgeMaterial));
      } else {
        // the piece only touches the hyperplane: its touching points
        const pts = new T3.BufferGeometry().setFromPoints(s.points.map((x) => new T3.Vector3(x[0], x[1], x[2])));
        add(new T3.Points(pts, new T3.PointsMaterial({ color: colour, size: 6, sizeAttenuation: false })));
      }
    });
  }

  function sayNoWebGL() {
    cv.title = "This drawing needs WebGL, which this browser does not offer.";
    const ctx = cv.getContext("2d");
    if (!ctx) return;
    const w = Math.max(1, cv.clientWidth);
    const hh = Math.max(1, cv.clientHeight);
    cv.width = w;
    cv.height = hh;
    ctx.fillStyle = cssVar("--muted") || "#666";
    ctx.font = "14px system-ui, sans-serif";
    ctx.textAlign = "center";
    ctx.fillText("This drawing needs WebGL, which this browser does not offer.", w / 2, hh / 2, w - 16);
  }

  function render() {
    if (destroyed || lost) return;
    if (failed) return sayNoWebGL();
    if (!THREE || !gl || !scene || !camera || !root || !last) return;
    build();
    const w = Math.max(1, cv.clientWidth);
    const hh = Math.max(1, cv.clientHeight);
    gl.setPixelRatio(window.devicePixelRatio || 1);
    gl.setSize(w, hh, false);
    // the whole container in view, whatever the slice and the turn: a sphere of the 4D container's circumradius (a 3D cut is inside it)
    const c = last.container;
    let r = c.kind === "ball" ? c.radius : 0;
    if (c.kind === "polytope") for (let v = 0; v < c.V.length / c.d; v++) r = Math.max(r, Math.hypot(...Array.from(c.V.subarray(v * c.d, (v + 1) * c.d))));
    const fov = (camera.fov * Math.PI) / 180;
    const fit = Math.min(1, w / hh);
    const dist = (1.12 * r) / Math.sin(Math.atan(Math.tan(fov / 2) * fit));
    camera.aspect = w / hh;
    camera.near = Math.max(1e-3, dist - 2 * r);
    camera.far = dist + 2 * r;
    camera.position.set(0, 0, dist);
    camera.lookAt(0, 0, 0);
    camera.updateProjectionMatrix();
    // about y by the yaw, then about x by the pitch, both in the view's frame
    root.rotation.set(pitch, yaw, 0, "XYZ");
    gl.render(scene, camera);
  }

  const onLost = () => (lost = true);
  const onRestored = () => {
    lost = false;
    render();
  };
  cv.addEventListener("webglcontextlost", onLost);
  cv.addEventListener("webglcontextrestored", onRestored);

  void Promise.all([import("three"), import("three/examples/jsm/geometries/ConvexGeometry.js")])
    .then(([T3, G]) => {
      if (destroyed) return;
      THREE = T3;
      Hull = G.ConvexGeometry;
      setup(T3);
      render();
    })
    .catch(() => {
      failed = true;
      render();
    });

  return {
    element: cv,
    attach(box: HTMLElement) {
      off = onDrag(cv, (dx, dy) => this.turn(dx * 0.01, dy * 0.01));
      void box;
    },
    draw(container: Placed, pieces: (Placed | null)[], axis: number, at: number, s: number) {
      last = { container, pieces, axis, at, s };
      render();
    },
    turn(dx: number, dy: number) {
      yaw += dx;
      pitch = Math.max(-Math.PI / 2, Math.min(Math.PI / 2, pitch + dy));
      render();
    },
    reset() {
      yaw = START.yaw;
      pitch = START.pitch;
      render();
    },
    destroy() {
      destroyed = true;
      off?.();
      cv.removeEventListener("webglcontextlost", onLost);
      cv.removeEventListener("webglcontextrestored", onRestored);
      dispose();
      if (gl) {
        gl.dispose();
        gl.forceContextLoss();
      }
      gl = null;
      scene = null;
    },
  };
}
