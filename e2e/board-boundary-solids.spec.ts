import { readFileSync } from "node:fs";
import { type Page, expect, test } from "@playwright/test";
import { type SpaceProblemResponse, validateSpaceSubmission } from "@polypack/spaces";
import { LAZY_VIEW, distinctColours, pageErrors } from "./helpers";

// The Platonic solids on the board (docs/spaces/BOUNDARY.md 4C): the picker offers the tetrahedron, octahedron, dodecahedron and icosahedron on
// every 4-polytope (the cube too, rotated except on the tesseract), the three.js views draw them, and a short run of the boundary search finds a
// packing that the server's own check accepts and that can be submitted (these problems have no research record: the first valid entry is kept).

declare const polypack: {
  space?: { lastResult: { config: unknown; score: number; valid: boolean; submittable: boolean } | null };
};

const SOLIDS = ["Balls", "Cubes", "Tetrahedra", "Octahedra", "Dodecahedra", "Icosahedra"];

/** The test packings of python/boundary/fixtures/platonic_packings.json, each certified strictly valid by the exact certifier. */
const PACKINGS = JSON.parse(readFileSync(new URL("../python/boundary/fixtures/platonic_packings.json", import.meta.url), "utf8")).packings as {
  id: string;
  config: { s: number; pieces: unknown[] };
}[];

/** Puts a test packing on the server as the first entry of its problem (a problem without a research record keeps the first valid one). */
async function seed(page: Page, id: string): Promise<void> {
  const p = PACKINGS.find((x) => x.id === id)!;
  const r = await page.request.post("/api/spaces/submit", { data: { space: "boundary", id, name: "E2E seed", units: "diameter", config: p.config } });
  expect(r.status(), await r.text()).toBe(200);
}

/** Waits until the board shows a 3D on 4D problem and its drawing is on screen. */
async function drawn(page: Page, id: string): Promise<void> {
  await expect(page.locator("#space-problem-note")).toHaveText(`Problem ${id}`, LAZY_VIEW);
  await expect(page.locator(".pp3 canvas")).toBeVisible(LAZY_VIEW);
}

/** Waits until the board shows a 3D on 4D problem that has no packing yet (the drawing says so). */
async function empty(page: Page, id: string): Promise<void> {
  await expect(page.locator("#space-problem-note")).toHaveText(`Problem ${id}`, LAZY_VIEW);
  await expect(page.locator(".pp3-status")).toHaveText("No packing to draw.", LAZY_VIEW);
}

const canvas = (page: Page) => page.locator(".pp3 canvas");
const button = (page: Page, name: string) => page.locator(".pp3").getByRole("button", { name, exact: true });

async function picture(page: Page): Promise<Buffer> {
  await page.evaluate(() => new Promise((r) => requestAnimationFrame(() => requestAnimationFrame(r))));
  return canvas(page).screenshot();
}

test("the picker offers the solids on every 4-polytope and the d-cubes only balls and cubes", async ({ page }) => {
  const errors = pageErrors(page);
  await page.goto("/#space=boundary/24cell/octahedron/3");
  await empty(page, "boundary/24cell/octahedron/3");
  await expect(page.locator("main")).toHaveAttribute("data-dimension", "3-on-4");
  await expect(page.locator("#space-f-piece")).toHaveValue("octahedron");
  await expect(page.locator("#space-f-piece option")).toHaveText(SOLIDS);
  await expect(page.locator("#space-f-P option")).toHaveCount(24);
  // no record: the lower bound is the problem's, the packing is none
  const api = (await (await page.request.get("/api/spaces/boundary/problems/boundary/24cell/octahedron/3")).json()) as SpaceProblemResponse;
  expect(api.status).toBe("no_record");
  expect(api.bound!.label).toMatch(/^\[CALC\] (volume|one piece)/);
  await expect(page.locator("#space-facts")).toContainText(`Lower bound${api.bound!.value.toFixed(6)}`);
  await expect(page.locator("#space-facts")).toContainText("no record yet");
  // every container takes every solid
  for (const polytope of ["tesseract", "5cell", "16cell", "120cell", "600cell"]) {
    await page.selectOption("#space-f-polytope", polytope);
    await expect(page.locator("#space-f-piece option")).toHaveText(SOLIDS);
    await expect(page.locator("#space-f-piece")).toHaveValue("octahedron");
  }
  await page.selectOption("#space-f-polytope", "600cell");
  await page.selectOption("#space-f-piece", "icosahedron");
  await empty(page, "boundary/600cell/icosahedron/3");
  await expect(page).toHaveURL(/#space=boundary\/600cell\/icosahedron\/3$/);
  expect(errors).toEqual([]);
});

test("a short run packs a solid, the three views draw it, and the server takes it as the first entry", async ({ page }) => {
  test.setTimeout(240_000);
  const errors = pageErrors(page);
  const id = "boundary/24cell/octahedron/1";
  await page.goto(`/#space=${id}`);
  await empty(page, id);
  await expect(page.locator("#space-method")).toContainText("Penalty relaxation of the pieces across the cells");
  await page.selectOption("#space-solve-budget", "evals");
  await page.fill("#space-solve-evals", "6000");
  await page.fill("#space-solve-seed", "1");
  await page.click("#space-solve-run");
  await expect(page.locator("#space-solve-stop")).toBeDisabled({ timeout: 150_000 });
  await expect(page.locator("#space-solve-run")).toBeEnabled();
  const result = (await page.evaluate(() => polypack.space!.lastResult))!;
  expect(result.valid).toBe(true);
  expect(result.submittable).toBe(true);
  // the server's own rule accepts it at its stated size, with a rotation the board kept
  const api = (await (await page.request.get(`/api/spaces/boundary/problems/${id}`)).json()) as SpaceProblemResponse;
  const v = validateSpaceSubmission("boundary", id, result.config, { bound: api.bound?.value ?? null, certified: null });
  expect(v.valid, v.reasons.join("; ")).toBe(true);
  expect(result.score).toBeGreaterThanOrEqual(api.bound!.value * (1 - 1e-9));
  await expect(page.locator(".pp3 canvas")).toBeVisible(LAZY_VIEW);
  await expect(page.locator(".pp3-status")).toContainText(/^1 octahedron, edge \d\.\d{6}\./);

  // drawn in Cells, as a net and in a slice: pieces on screen each time
  await expect(button(page, "Cells")).toHaveAttribute("aria-pressed", "true");
  const cells = await picture(page);
  expect(await distinctColours(page, cells)).toBeGreaterThan(30);
  await button(page, "Net").click();
  await expect(page.getByRole("slider", { name: "Fold" })).toBeVisible();
  const net = await picture(page);
  expect(net.equals(cells)).toBe(false);
  expect(await distinctColours(page, net)).toBeGreaterThan(30);
  await button(page, "Slice").click();
  await expect(page.getByRole("slider", { name: "Position" })).toBeVisible();
  expect(await distinctColours(page, await picture(page))).toBeGreaterThan(20);

  // the first valid entry of a problem without a record is kept
  await expect(page.locator("#space-submit")).toBeEnabled();
  await page.fill("#name", "E2E solids");
  const sent = page.waitForRequest((r) => r.url().endsWith("/api/spaces/submit") && r.method() === "POST");
  const answered = page.waitForResponse((r) => r.url().endsWith("/api/spaces/submit"));
  await page.click("#space-submit");
  const body = (await sent).postDataJSON();
  expect(body).toMatchObject({ space: "boundary", id, name: "E2E solids" });
  expect(body.config.pieces[0].q).toHaveLength(4);
  expect((await answered).status()).toBe(200);
  await expect(page.locator("#space-submit-result")).toContainText("✓ Kept. Server score");
  expect(errors).toEqual([]);
});

test("another solid is drawn in the view last used, and the dimensions page counts the solids", async ({ page }) => {
  const errors = pageErrors(page);
  // two test packings, as site entries (the server's rule judges them first)
  await seed(page, "boundary/tesseract/icosahedron/2");
  await seed(page, "boundary/tesseract/tetrahedron/3");
  await page.goto("/#space=boundary/tesseract/icosahedron/2");
  await drawn(page, "boundary/tesseract/icosahedron/2");
  await expect(page.locator(".pp3-status")).toContainText(/^2 icosahedra, edge 2\.1376/);
  await button(page, "Net").click();
  await expect(button(page, "Net")).toHaveAttribute("aria-pressed", "true");
  await page.selectOption("#space-f-piece", "tetrahedron");
  await page.selectOption("#space-f-P", "3");
  await drawn(page, "boundary/tesseract/tetrahedron/3");
  await expect(page.locator(".pp3-status")).toContainText(/^3 tetrahedra, edge 0\.5772/);
  await expect(button(page, "Net")).toHaveAttribute("aria-pressed", "true");
  expect(await distinctColours(page, await picture(page))).toBeGreaterThan(30);
  // the dimensions page counts the 3D on 4D problems and opens the table of the solids
  await page.goto("/dimensions.html");
  await expect(page.locator('#dimension-groups > section[data-dimension="3-on-4"]')).toContainText("816");
  expect(errors).toEqual([]);
});
