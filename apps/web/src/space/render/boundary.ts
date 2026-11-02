/**
 * The board's drawing of a 3D on 4D problem: the three.js views of packages/view3d (Cells, Net
 * with its Fold slider, Slice with its axis and position), wrapped in the board's SpaceView.
 *
 * three.js is large, so it is loaded with a dynamic import when the first boundary problem is
 * drawn, not with the board's space mode: the other kinds never load it. Until it arrives, the
 * drawing says that it is loading; a configuration shown meanwhile is drawn when it is ready.
 *
 * One view serves every problem of the same polytope and piece (canShow), so that choosing
 * another P keeps the view, the camera, the fold and the slice the user chose. Another piece
 * needs another view (the piece's shape is fixed in a view); it opens in the view last used.
 */
import type { BoundaryView, ViewKey } from "@polypack/view3d";
import type { BoundaryConfig, PieceName, PolytopeName } from "@polypack/boundary";
import type { SpaceConfig } from "@polypack/spaces";
import { pieceNames } from "../../boundary-model";
import type { SpaceView } from "../view";

/** The fields of a boundary problem the drawing needs. */
export interface BoundaryDrawing {
  polytope: string;
  piece: string;
}

type View3d = typeof import("@polypack/view3d");

let module: Promise<View3d> | null = null;

/** Loads packages/view3d (and three.js) once; a failed load is tried again on the next call. */
export function loadView3d(): Promise<View3d> {
  return (module ??= import("@polypack/view3d").catch((error: unknown) => {
    module = null;
    throw error;
  }));
}

/** The view last shown (Cells, Net or Slice), for the next drawing of another piece. */
let lastView: ViewKey = "atlas";

export function mountBoundaryView(host: HTMLElement, p: BoundaryDrawing, label = "Drawing"): SpaceView {
  const root = document.createElement("div");
  root.className = "space-view space-view-3d";
  const wait = document.createElement("div");
  wait.className = "space-3d-wait";
  const text = document.createElement("p");
  text.setAttribute("role", "status");
  text.textContent = "Loading the 3D drawing…";
  wait.append(text);
  root.append(wait);
  host.replaceChildren(root);

  // The d-cubes (5cube to 10cube, 4D on 5D to 9D on 10D) are drawn as a table and a chart (ncube.ts), and so are the simplex, orthoplex and
  // demicube of the same dimensions (ndpoly.ts), not as the 3D views of the tesseract.
  if (/^\d+(cube|simplex|orthoplex|demicube)$/.test(p.polytope)) {
    let inner: SpaceView | null = null;
    let gone = false;
    let later: SpaceConfig | null | undefined;
    const load = /^\d+cube$/.test(p.polytope)
      ? import("./ncube").then((m) => m.mountNcubeView(host, p as Parameters<typeof m.mountNcubeView>[1], label))
      : import("./ndpoly").then((m) => m.mountNdpolyView(host, p as Parameters<typeof m.mountNdpolyView>[1], label));
    void load.then(
      (view) => {
        if (gone) {
          view.destroy();
          return;
        }
        inner = view;
        if (later !== undefined) inner.show(later);
      },
      () => {
        text.textContent = "The drawing could not be loaded. The connection may have dropped.";
      },
    );
    return {
      show(c) {
        if (inner) inner.show(c);
        else later = c;
      },
      canShow: (q) => q.space === "boundary" && q.polytope === p.polytope && q.piece === p.piece,
      destroy() {
        gone = true;
        inner?.destroy();
        root.remove();
      },
    };
  }

  let view: BoundaryView | null = null;
  let gone = false;
  /** The configuration to draw once the view exists (undefined: none was asked for). */
  let pending: SpaceConfig | null | undefined;

  function open(): void {
    text.textContent = "Loading the 3D drawing…";
    wait.querySelector("button")?.remove();
    loadView3d().then(
      (m) => {
        if (gone) return;
        root.replaceChildren();
        view = m.createBoundaryView(root, {
          ...m.engineFor(p.polytope as PolytopeName, p.piece as PieceName),
          initialView: lastView,
          label,
          pieceNames: pieceNames(p.piece),
        });
        if (pending !== undefined) view.show(pending as BoundaryConfig | null);
      },
      () => {
        if (gone) return;
        text.textContent = "The 3D drawing could not be loaded. The connection may have dropped.";
        const retry = document.createElement("button");
        retry.type = "button";
        retry.textContent = "Try again";
        retry.addEventListener("click", open);
        wait.append(retry);
      },
    );
  }
  open();

  return {
    show(config) {
      // a configuration of this problem's format (read by the space's reader, model.ts)
      if (view) view.show(config as BoundaryConfig | null);
      else pending = config;
    },
    canShow(q) {
      return q.space === "boundary" && q.polytope === p.polytope && q.piece === p.piece;
    },
    destroy() {
      gone = true;
      if (view) {
        lastView = view.view;
        view.destroy();
        view = null;
      }
      root.remove();
    },
  };
}
