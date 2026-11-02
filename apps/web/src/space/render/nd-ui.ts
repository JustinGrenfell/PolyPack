/**
 * Small DOM helpers shared by the board's drawings of the nD spaces (ndim.ts: 4D in 4D to 10D in 10D; ncube.ts: 4D on 5D to
 * 9D on 10D): an element builder, a labelled select, a table of pieces and an SVG chart frame. Text is always inserted
 * as text.
 */
import { SVG_NS, cssVar, isDark } from "./common";

export function h<K extends keyof HTMLElementTagNameMap>(
  tag: K,
  attrs: Record<string, string | boolean> = {},
  ...children: (Node | string)[]
): HTMLElementTagNameMap[K] {
  const e = document.createElement(tag);
  for (const [k, v] of Object.entries(attrs)) {
    if (v === false) continue;
    if (k === "class") e.className = String(v);
    else e.setAttribute(k, v === true ? "" : v);
  }
  e.append(...children);
  return e;
}

export function svg<K extends keyof SVGElementTagNameMap>(
  tag: K,
  attrs: Record<string, string | number> = {},
  ...children: (Node | string)[]
): SVGElementTagNameMap[K] {
  const e = document.createElementNS(SVG_NS, tag);
  for (const [k, v] of Object.entries(attrs)) e.setAttribute(k, String(v));
  e.append(...children);
  return e;
}

let counter = 0;

/** A select with a visible label, in a toolbar. `options` are [value, text]. */
export function labelledSelect(
  label: string,
  options: readonly [string, string][],
  value: string,
  onChange: (value: string) => void,
  attrs: Record<string, string> = {},
): { box: HTMLElement; select: HTMLSelectElement } {
  const id = `nd-select-${++counter}`;
  const select = h("select", { id, ...attrs });
  for (const [v, t] of options) select.append(new Option(t, v));
  select.value = value;
  select.addEventListener("change", () => onChange(select.value));
  const box = h("span", { class: "nd-field" }, h("label", { for: id }, label), select);
  return { box, select };
}

/** A fill for piece i, in the page's hues (common.ts pieceFill), for SVG. */
export { pieceFill } from "./common";

/** The page's colours for a chart: the container's edge and face, and the text. */
export function chartColours(): { edge: string; face: string; text: string; piece: string } {
  const dark = isDark();
  return {
    edge: cssVar("--space-edge") || (dark ? "#c9cbd0" : "#3b3b3f"),
    face: cssVar("--space-face") || (dark ? "#2a2c31" : "#ece8df"),
    text: cssVar("--muted") || (dark ? "#aaa" : "#555"),
    piece: cssVar("--space-piece-edge") || (dark ? "rgba(255,255,255,0.5)" : "rgba(0,0,0,0.45)"),
  };
}

/**
 * A table of rows under a caption, in a box that scrolls sideways on a narrow screen. `head` and every row are text
 * cells; a cell may also be a node. Numeric columns (by index) are right-aligned by the `num` class.
 */
export function dataTable(
  caption: string,
  head: readonly string[],
  rows: readonly (readonly (string | Node)[])[],
  numeric: readonly number[] = [],
): HTMLElement {
  const table = h("table", { class: "nd-table" });
  table.append(h("caption", {}, caption));
  table.append(h("thead", {}, h("tr", {}, ...head.map((t) => h("th", { scope: "col" }, t)))));
  const body = h("tbody");
  for (const r of rows)
    body.append(
      h(
        "tr",
        {},
        ...r.map((c, i) => h("td", numeric.includes(i) ? { class: "num" } : {}, typeof c === "string" ? c : c)),
      ),
    );
  table.append(body);
  return h("div", { class: "nd-table-wrap", tabindex: "0", role: "region", "aria-label": caption }, table);
}
