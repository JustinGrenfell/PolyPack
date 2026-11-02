/**
 * The board's drawing of a (d-1)D on dD packing on the simplex, orthoplex or demicube, d = 5 to 10 (4D on 5D to 9D on 10D):
 * balls on the boundary of the container (boundary/<d><simplex|orthoplex|demicube>/ball/<P>, docs/spaces/BOUNDARY.md 4C). A
 * configuration is { s, pieces: [{ cell, u }] }: the cell is a facet of the container (named by `facetLabel`) and u the ball's centre in
 * the facet's frame, its centre at the origin.
 *
 * What is drawn is light, as the d-cubes' (ncube.ts): a table of the balls (index, facet, centre in the facet, centre in the engine's
 * coordinates x1 to xD) and a chart of the centres, placed in R^D and projected onto two axes the user picks, in the colour of their
 * facet, inside the projected container (the outline of its projected vertices). A ball on a demicube facet x_i = const is drawn as a disk
 * when neither axis is i; the others as a dot.
 */
import { type NDName, getND, ndParse } from "@polypack/boundary";
import type { SpaceConfig } from "@polypack/spaces";
import type { SearchProblem } from "@polypack/spaces/search";
import type { SpaceView } from "../view";
import { pieceFill } from "./common";
import { centreText } from "./ndslice";
import { chartColours, dataTable, h, labelledSelect, svg } from "./nd-ui";

type BoundaryProblem = Extract<SearchProblem, { space: "boundary" }>;

interface NPiece {
  cell: number;
  u: number[];
}

/** The name of a container of this view ("7simplex"), or null for the others (the d-cubes, the 4-polytopes). */
export const ndpolyName = (polytope: string): NDName | null => (ndParse(polytope) ? (polytope as NDName) : null);

/** A configuration read for the drawing, or null when it is not one of this container: right piece form, finite numbers. */
export function readNdpolyConfig(c: unknown, polytope: NDName): { s: number; pieces: NPiece[] } | null {
  if (typeof c !== "object" || c === null) return null;
  const { s, pieces } = c as { s?: unknown; pieces?: unknown };
  const C = getND(polytope);
  if (typeof s !== "number" || !Number.isFinite(s) || s <= 0 || !Array.isArray(pieces)) return null;
  const out: NPiece[] = [];
  for (const q of pieces) {
    if (typeof q !== "object" || q === null) return null;
    const { cell, u } = q as { cell?: unknown; u?: unknown };
    if (!Number.isInteger(cell) || (cell as number) < 0 || (cell as number) >= C.C) return null;
    if (!Array.isArray(u) || u.length !== C.n || !u.every((x) => typeof x === "number" && Number.isFinite(x))) return null;
    out.push({ cell: cell as number, u: u as number[] });
  }
  return { s, pieces: out };
}

/** The piece's centre in the engine's coordinates R^D. */
export function ndWorldCentre(polytope: NDName, q: NPiece, s: number): number[] {
  return Array.from(getND(polytope).toAmbient(q.cell, q.u, s));
}

/** The convex hull of 2D points (monotone chain), counter-clockwise. */
export function hull2(points: readonly [number, number][]): [number, number][] {
  const p = [...points].sort((a, b) => a[0] - b[0] || a[1] - b[1]);
  if (p.length < 3) return p;
  const cross = (o: [number, number], a: [number, number], b: [number, number]) => (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0]);
  const lo: [number, number][] = [];
  for (const q of p) {
    while (lo.length >= 2 && cross(lo[lo.length - 2], lo[lo.length - 1], q) <= 0) lo.pop();
    lo.push(q);
  }
  const up: [number, number][] = [];
  for (const q of [...p].reverse()) {
    while (up.length >= 2 && cross(up[up.length - 2], up[up.length - 1], q) <= 0) up.pop();
    up.push(q);
  }
  return [...lo.slice(0, -1), ...up.slice(0, -1)];
}

export function mountNdpolyView(host: HTMLElement, p: BoundaryProblem, label = "Drawing"): SpaceView {
  const name = ndpolyName(p.polytope)!;
  const C = getND(name);
  const D = C.D;
  const info = ndParse(name)!;
  const root = h("div", { class: "space-view space-view-nd space-view-nd-chart" });
  const tools = h("div", { class: "toolbar space-tools nd-tools", role: "toolbar", "aria-label": label });
  const status = h("p", { class: "note nd-status", role: "status" }, "No packing to draw.");
  const box = h("div", { class: "space-draw nd-draw", role: "img" });
  const tableHost = h("div", { class: "nd-pieces" });
  const summary = h("p", { class: "note nd-facets" });
  const element = svg("svg", { class: "nd-chart", role: "presentation", preserveAspectRatio: "xMidYMid meet" });
  box.append(element);
  root.append(tools, status, box, summary, tableHost);
  host.replaceChildren(root);

  const axes = Array.from({ length: D }, (_, k): [string, string] => [String(k), `x${k + 1}`]);
  let ax = 0;
  let ay = 1;
  let config: { s: number; pieces: NPiece[] } | null = null;
  const fx = labelledSelect("Horizontal axis", axes, String(ax), (v) => {
    ax = Number(v);
    if (ax === ay) {
      ay = (ax + 1) % D;
      fy.select.value = String(ay);
    }
    draw();
  }, { "data-nd": "axis-x" });
  const fy = labelledSelect("Vertical axis", axes, String(ay), (v) => {
    ay = Number(v);
    if (ax === ay) {
      ax = (ay + 1) % D;
      fx.select.value = String(ax);
    }
    draw();
  }, { "data-nd": "axis-y" });
  tools.append(fx.box, fy.box);

  function draw(): void {
    if (!config) {
      status.hidden = false;
      box.style.visibility = "hidden";
      box.setAttribute("aria-label", "No packing to draw.");
      return;
    }
    status.hidden = true;
    box.style.visibility = "";
    const { s, pieces } = config;
    const col = chartColours();
    const outline = hull2(C.vertices(s).map((v): [number, number] => [v[ax], -v[ay]]));
    const centres = pieces.map((q) => ndWorldCentre(name, q, s));
    const xs = [...outline.map((q) => q[0]), ...centres.map((c) => c[ax] - 0.5), ...centres.map((c) => c[ax] + 0.5)];
    const ys = [...outline.map((q) => q[1]), ...centres.map((c) => -c[ay] - 0.5), ...centres.map((c) => -c[ay] + 0.5)];
    const x0 = Math.min(...xs);
    const x1 = Math.max(...xs);
    const y0 = Math.min(...ys);
    const y1 = Math.max(...ys);
    const span = Math.max(x1 - x0, y1 - y0, 1e-9);
    const m = 0.08 * span;
    element.setAttribute("viewBox", `${x0 - m} ${y0 - m} ${x1 - x0 + 2 * m} ${y1 - y0 + 2 * m}`);
    const stroke = { "vector-effect": "non-scaling-stroke", "stroke-linejoin": "round" } as const;
    const parts: SVGElement[] = [
      svg("polygon", { points: outline.map((q) => q.join(",")).join(" "), fill: col.face, "fill-opacity": 0.5, stroke: col.edge, "stroke-width": 1.5, ...stroke }),
    ];
    pieces.forEach((q, n) => {
      const c = centres[n];
      const fill = pieceFill(q.cell, 0.55);
      const g = svg("g", { "data-piece": n + 1, "data-cell": q.cell });
      // a ball in a demicube facet x_i = const projects to a full disk on two other axes
      const full = info.family === "demicube" && q.cell < 2 * info.d && Math.floor(q.cell / 2) !== ax && Math.floor(q.cell / 2) !== ay;
      if (full) g.append(svg("circle", { cx: c[ax], cy: -c[ay], r: 0.5, fill, stroke: col.piece, "stroke-width": 1, ...stroke }));
      else g.append(svg("circle", { cx: c[ax], cy: -c[ay], r: 0.04 * span, fill, stroke: col.piece, "stroke-width": 1, ...stroke }));
      g.append(
        svg("text", { x: c[ax], y: -c[ay], "font-size": Math.max(0.025 * span, full ? 0.3 : 0), "text-anchor": "middle", "dominant-baseline": "central", fill: col.text }, String(n + 1)),
      );
      parts.push(g);
    });
    element.replaceChildren(...parts);
    box.setAttribute(
      "aria-label",
      `${pieces.length} balls on the boundary of the ${info.d}-${info.family} of edge ${s.toFixed(6)}, their centres projected onto x${ax + 1} and x${ay + 1}.`,
    );
  }

  function renderTable(): void {
    if (!config) {
      tableHost.replaceChildren();
      summary.textContent = "";
      return;
    }
    const counts = new Map<number, number>();
    for (const q of config.pieces) counts.set(q.cell, (counts.get(q.cell) ?? 0) + 1);
    const used = [...counts.entries()].sort((a, b) => a[0] - b[0]);
    summary.textContent = `Balls per facet: ${used.map(([c, n]) => `${C.facetLabel(c)} ${n}`).join(", ") || "none"}. The other ${C.C - used.length} of the ${C.C} facets hold none.`;
    const rows = config.pieces.map((q, i) => [String(i + 1), C.facetLabel(q.cell), centreText(q.u), centreText(ndWorldCentre(name, q, config!.s))]);
    tableHost.replaceChildren(
      dataTable(
        `Balls of the packing (${config.pieces.length}); a centre in its facet is in the facet's own frame, in the engine's coordinates x1 to x${D}`,
        ["#", "Facet", "Centre in the facet", `Centre in x1 to x${D}`],
        rows,
        [0],
      ),
    );
  }

  const resize = new ResizeObserver(() => draw());
  resize.observe(box);
  const scheme = matchMedia("(prefers-color-scheme: dark)");
  scheme.addEventListener("change", draw);

  return {
    show(c: SpaceConfig | null) {
      config = c ? readNdpolyConfig(c, name) : null;
      renderTable();
      draw();
    },
    canShow: (q) => q.space === "boundary" && q.polytope === p.polytope && q.piece === p.piece,
    destroy() {
      resize.disconnect();
      scheme.removeEventListener("change", draw);
      root.remove();
    },
  };
}
