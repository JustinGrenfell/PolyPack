/**
 * The board's drawing of a (d-1)D on dD packing, 4D on 5D to 9D on 10D: balls or unit cubes on the boundary of the d-cube
 * (boundary/<d>cube/<ball|cube>/<P>, docs/spaces/BOUNDARY.md 4B). A configuration is { s, pieces: [{ cell, u }] }: cell
 * 2i + (sigma > 0) is the facet x_i = sigma s/2 (a (d-1)-cube), and u the piece's centre in that facet's coordinates, the
 * other d-1 coordinates in increasing order.
 *
 * What is drawn is light: a table of the pieces (index, facet, centre in the facet) and a chart of the centres, placed in R^d
 * on the facet and projected onto two axes the user picks, in the colour of their facet, inside the projected d-cube (a
 * square). A ball is drawn as a disk when both axes lie in its facet.
 */
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

/** "d" of a container name "7cube". */
export const ncubeDimension = (polytope: string): number | null => {
  const m = /^(\d+)cube$/.exec(polytope);
  return m ? Number(m[1]) : null;
};

/** "+x3" or "−x1": the facet of a cell, 1-based axes. */
export const facetName = (cell: number): string => `${cell % 2 ? "+" : "−"}x${Math.floor(cell / 2) + 1}`;

/** A configuration read for the drawing, or null when it is not one of this container: right piece form, finite numbers. */
export function readNcubeConfig(c: unknown, d: number): { s: number; pieces: NPiece[] } | null {
  if (typeof c !== "object" || c === null) return null;
  const { s, pieces } = c as { s?: unknown; pieces?: unknown };
  if (typeof s !== "number" || !Number.isFinite(s) || s <= 0 || !Array.isArray(pieces)) return null;
  const out: NPiece[] = [];
  for (const q of pieces) {
    if (typeof q !== "object" || q === null) return null;
    const { cell, u } = q as { cell?: unknown; u?: unknown };
    if (!Number.isInteger(cell) || (cell as number) < 0 || (cell as number) >= 2 * d) return null;
    if (!Array.isArray(u) || u.length !== d - 1 || !u.every((x) => typeof x === "number" && Number.isFinite(x))) return null;
    out.push({ cell: cell as number, u: u as number[] });
  }
  return { s, pieces: out };
}

/** The piece's centre in R^d: the facet's axis at +-s/2 and the other coordinates from u. */
export function worldCentre(q: NPiece, d: number, s: number): number[] {
  const axis = Math.floor(q.cell / 2);
  const out = new Array<number>(d);
  let k = 0;
  for (let i = 0; i < d; i++) out[i] = i === axis ? (q.cell % 2 ? s / 2 : -s / 2) : q.u[k++];
  return out;
}

export function mountNcubeView(host: HTMLElement, p: BoundaryProblem, label = "Drawing"): SpaceView {
  const d = ncubeDimension(p.polytope)!;
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

  const axes = Array.from({ length: d }, (_, k): [string, string] => [String(k), `x${k + 1}`]);
  let ax = 0;
  let ay = 1;
  let config: { s: number; pieces: NPiece[] } | null = null;
  // until the user picks axes, the chart uses the two axes that the fewest pieces have their facet on, so that the pieces
  // of the other facets show in full instead of as dots on the edge
  let touched = false;
  const fx = labelledSelect("Horizontal axis", axes, String(ax), (v) => {
    touched = true;
    ax = Number(v);
    if (ax === ay) {
      ay = (ax + 1) % d;
      fy.select.value = String(ay);
    }
    draw();
  }, { "data-nd": "axis-x" });
  const fy = labelledSelect("Vertical axis", axes, String(ay), (v) => {
    touched = true;
    ay = Number(v);
    if (ax === ay) {
      ax = (ay + 1) % d;
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
    const hs = s / 2;
    const col = chartColours();
    const span = s;
    const m = 0.1 * span;
    element.setAttribute("viewBox", `${-hs - m} ${-hs - m} ${s + 2 * m} ${s + 2 * m}`);
    const stroke = { "vector-effect": "non-scaling-stroke", "stroke-linejoin": "round" } as const;
    const parts: SVGElement[] = [
      svg("rect", { x: -hs, y: -hs, width: s, height: s, fill: col.face, "fill-opacity": 0.5, stroke: col.edge, "stroke-width": 1.5, ...stroke }),
    ];
    pieces.forEach((q, n) => {
      const c = worldCentre(q, d, s);
      const facetAxis = Math.floor(q.cell / 2);
      const inFacet = facetAxis !== ax && facetAxis !== ay;
      // the facet's hue: the same for its pieces, so that the pieces of a facet read as a group
      const fill = pieceFill(q.cell, 0.55);
      const g = svg("g", { "data-piece": n + 1, "data-cell": q.cell });
      if (p.piece === "ball" && inFacet)
        g.append(svg("circle", { cx: c[ax], cy: -c[ay], r: 0.5, fill, stroke: col.piece, "stroke-width": 1, ...stroke }));
      else if (p.piece === "cube" && inFacet)
        g.append(svg("rect", { x: c[ax] - 0.5, y: -c[ay] - 0.5, width: 1, height: 1, fill, stroke: col.piece, "stroke-width": 1, ...stroke }));
      else g.append(svg("circle", { cx: c[ax], cy: -c[ay], r: 0.08 * Math.min(1, span), fill, stroke: col.piece, "stroke-width": 1, ...stroke }));
      g.append(
        svg(
          "text",
          { x: c[ax], y: -c[ay], "font-size": 0.035 * span, "text-anchor": "middle", "dominant-baseline": "central", fill: col.text },
          String(n + 1),
        ),
      );
      parts.push(g);
    });
    element.replaceChildren(...parts);
    box.setAttribute(
      "aria-label",
      `${pieces.length} ${p.piece === "ball" ? "balls" : "cubes"} on the boundary of the ${d}-cube of edge ${s.toFixed(6)}, their centres projected onto x${ax + 1} and x${ay + 1}.`,
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
    summary.textContent = `Pieces per facet: ${used.map(([c, n]) => `${facetName(c)} ${n}`).join(", ") || "none"}. The other ${2 * d - used.length} of the ${2 * d} facets hold none.`;
    const rows = config.pieces.map((q, i) => [String(i + 1), facetName(q.cell), centreText(q.u), centreText(worldCentre(q, d, config!.s))]);
    tableHost.replaceChildren(
      dataTable(
        `Pieces of the packing (${config.pieces.length}); a facet's centre is in its own coordinates, the world centre in x1 to x${d}`,
        ["#", "Facet", "Centre in the facet", `Centre in x1 to x${d}`],
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
      config = c ? readNcubeConfig(c, d) : null;
      if (config && !touched) {
        const use = Array.from({ length: d }, () => 0);
        for (const q of config.pieces) use[Math.floor(q.cell / 2)]++;
        const order = use.map((n, k) => [n, k]).sort((a, b) => a[0] - b[0] || a[1] - b[1]);
        [ax, ay] = [order[0][1], order[1][1]].sort((a, b) => a - b);
        fx.select.value = String(ax);
        fy.select.value = String(ay);
      }
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
