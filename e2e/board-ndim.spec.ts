import AxeBuilder from "@axe-core/playwright";
import { type Page, expect, test } from "@playwright/test";
import type { SpaceProblemResponse } from "@polypack/spaces";
import { DIMENSIONS, familyKey } from "../apps/web/src/dimension-table";
import { LAZY_VIEW, distinctColours, openProblem, pageErrors } from "./helpers";

// The nD dimensions on the board against the local stack (playwright.config.ts): 4D in 4D as a 3D slice, 5D in 5D to 10D in
// 10D as a chart and a table of the pieces, 4D on 5D to 9D on 10D as a chart of centres per facet and a table; the Records
// page of one dimension. The phone layout is in board-ndim-mobile.spec.ts.

/** Waits until the space mode shows a problem and its values. */
async function shown(page: Page, id: string, wait: { timeout?: number } = {}): Promise<void> {
  await expect(page.locator("#space-problem-note")).toHaveText(`Problem ${id}`, wait);
  await expect
    .poll(
      async () =>
        (await page.locator("#space-empty").isVisible()) ||
        ((await page.locator("#space-score").textContent()) ?? "—").trim() !== "—",
      wait,
    )
    .toBe(true);
}

/** Runs the optimizer with an evaluation budget and waits for its end; returns the status line. */
async function runFor(page: Page, evals: number, seed = 1): Promise<string> {
  await page.selectOption("#space-solve-budget", "evals");
  await page.fill("#space-solve-evals", String(evals));
  await page.fill("#space-solve-seed", String(seed));
  await page.click("#space-solve-run");
  await expect(page.locator("#space-solve-stop")).toBeDisabled({ timeout: 90_000 });
  await expect(page.locator("#space-solve-run")).toBeEnabled();
  return (await page.locator("#space-solve-status").textContent()) ?? "";
}

test("the Problem set menu lists the nD dimensions after 3D on 4D, and each opens its own dimension", async ({ page }) => {
  // modules of the nD drawings and their data, by their paths in the dev server or their file names in a build
  const nd = /(render\/(ndim|ncube|nd-ui|ndslice)|ConvexGeometry|data\/d[467]\.|\/(ndim|ncube|nd-ui|d4|d6|d7)-[A-Za-z0-9_-]{8}\.js)/;
  const loaded: string[] = [];
  page.on("request", (r) => {
    if (nd.test(new URL(r.url()).pathname)) loaded.push(r.url());
  });
  await openProblem(page, 4, 4, 5);
  const menu = page.locator("#family-select");
  await expect(menu.locator("optgroup")).toHaveCount(DIMENSIONS.length);
  const labels = await menu.locator("optgroup").evaluateAll((g) => g.map((x) => x.getAttribute("label")));
  expect(labels.slice(3, 8)).toEqual(["3D on 4D", "4D in 4D", "4D on 5D", "5D in 5D", "5D on 6D"]);
  expect(labels.at(-1)).toBe("10D in 10D");
  await expect(menu.locator('option[value="ndim-4"]')).toHaveText("4D in 4D: Polytopes and balls");
  await expect(menu.locator('option[value="boundary-5cube"]')).toHaveText("4D on 5D: 5-cube boundary");
  await expect(menu.locator("option[disabled]")).toHaveCount(0);

  // the default board has loaded none of it, and neither has another kind of the board
  expect(loaded).toEqual([]);
  await page.selectOption("#family-select", "sphere");
  await shown(page, "sphere/4/12", LAZY_VIEW);
  expect(loaded).toEqual([]);

  // a dimension's option opens its problem, in that dimension
  await page.selectOption("#family-select", "ndim-6");
  await shown(page, "ndim/6/cube/cube/12", LAZY_VIEW);
  await expect(page.locator("main")).toHaveAttribute("data-dimension", "6-in-6");
  await expect(page).toHaveURL(/#space=ndim\/6\/cube\/cube\/12$/);
  await expect(page.locator("#family-select option:checked")).toHaveText("6D in 6D: Polytopes and balls");
  // the picker offers 6D's shapes and range of P, and does not show the dimension
  await expect(page.locator("#space-f-d")).toHaveCount(0);
  await expect(page.locator("#space-f-piece option")).toHaveText([
    "Simplex",
    "Hypercube",
    "Orthoplex",
    "Demicube",
    "2_21 polytope",
    "1_22 polytope",
    "Ball",
  ]);
  await expect(page.locator("#space-f-P option")).toHaveCount(23);
  await expect(page.locator("#space-score-label")).toHaveText("Container scale s");
  expect(loaded.length).toBeGreaterThan(0);

  await page.selectOption("#family-select", "boundary-7cube");
  await shown(page, "boundary/7cube/ball/9", LAZY_VIEW);
  await expect(page.locator("main")).toHaveAttribute("data-dimension", "6-on-7");
  await expect(page.locator("#space-f-polytope")).toHaveCount(0);
  await expect(page.locator("#space-f-piece option")).toHaveText(["Balls", "Cubes"]);
  await expect(page.locator("#space-score-label")).toHaveText("7-cube edge s");

  // a family chosen again returns to the problem last shown in it; another one starts on its own
  await page.selectOption("#family-select", "ndim-6");
  await shown(page, "ndim/6/cube/cube/12");
  await page.selectOption("#space-f-container", "ortho");
  await shown(page, "ndim/6/cube/ortho/12");
  await page.selectOption("#family-select", "ndim-9");
  await shown(page, "ndim/9/cube/cube/12");
  await page.selectOption("#family-select", "ndim-6");
  await shown(page, "ndim/6/cube/ortho/12");

  // 3D on 4D keeps its option's value (the space), and the Flat board comes back as it was
  await page.selectOption("#family-select", "boundary");
  await shown(page, "boundary/tesseract/ball/9", LAZY_VIEW);
  await expect(page.locator("main")).toHaveAttribute("data-dimension", "3-on-4");
  await expect(page.locator("#space-f-polytope option")).toHaveCount(6);
  await page.selectOption("#family-select", "polygon");
  await expect(page.locator("#board")).toBeVisible();
});

test("addresses open the dimension that holds their problem, and the older ones still work", async ({ page }) => {
  const errors = pageErrors(page);
  for (const [hash, id, key, dimension] of [
    ["#space=ndim/7/ball/321/9", "ndim/7/ball/321/9", "ndim-7", "7-in-7"],
    ["#space=ndim", "ndim/4/cube/cube/12", "ndim-4", "4-in-4"],
    ["#space=boundary/10cube/cube/20", "boundary/10cube/cube/20", "boundary-10cube", "9-on-10"],
    ["#space=boundary/tesseract/cube/5", "boundary/tesseract/cube/5", "boundary", "3-on-4"],
    ["#space=boundary/600cell/ball/2", "boundary/600cell/ball/2", "boundary", "3-on-4"],
    ["#space=boundary", "boundary/tesseract/ball/9", "boundary", "3-on-4"],
  ] as const) {
    await page.goto("/");
    await page.goto(`/${hash}`);
    await page.reload();
    await shown(page, id, LAZY_VIEW);
    await expect(page.locator("#family-select")).toHaveValue(key);
    await expect(page.locator("main")).toHaveAttribute("data-dimension", dimension);
  }
  // a hash change inside the page does the same, and Back to the Flat board works
  await page.evaluate(() => (location.hash = "#space=ndim/5/ortho/ortho/4"));
  await shown(page, "ndim/5/ortho/ortho/4");
  await expect(page.locator("#family-select")).toHaveValue("ndim-5");
  await page.evaluate(() => (location.hash = "#p=4,4,5"));
  await expect(page.locator("#board")).toBeVisible();
  await expect(page.locator("#family-select")).toHaveValue("polygon");
  expect(errors).toEqual([]);
});

test("4D in 4D draws a 3D slice, cut by the axis and the position the user picks", async ({ page }) => {
  const errors = pageErrors(page);
  await page.goto("/");
  await page.goto("/#space=ndim/4/cube/ortho/12");
  await page.reload();
  await shown(page, "ndim/4/cube/ortho/12", LAZY_VIEW);
  const draw = page.locator(".nd-draw");
  await expect(draw.locator("canvas")).toBeVisible(LAZY_VIEW);
  await expect(draw).toHaveAttribute("aria-label", /^12 hypercubes in an orthoplex of scale 4\.000000, cut by x4 = 0\.000, drawn in 3D\.$/);
  // the research record, with its values from the API
  const api = (await (await page.request.get("/api/spaces/ndim/problems/ndim/4/cube/ortho/12")).json()) as SpaceProblemResponse;
  await expect(page.locator("#space-facts")).toContainText(`Certified value${api.certified!.display}`);
  await expect(page.locator("#space-shown")).toHaveText("the research record");
  // WebGL drew something besides the clear colour
  await expect.poll(async () => distinctColours(page, await draw.locator("canvas").screenshot())).toBeGreaterThan(8);

  // the slider covers the container's extent along the axis; the readout and the drawing follow it
  const slider = page.locator('[data-nd="position"]');
  const lo = Number(await slider.getAttribute("min"));
  const hi = Number(await slider.getAttribute("max"));
  // the 16-cell's vertices are +-e_i / sqrt 2, scaled by s = 4: the extent along an axis is +-2 sqrt 2
  // (the record's scale is 4 and a hair, 4.0000002, which keeps its 1e-7 clearance)
  expect(lo).toBeCloseTo(-2 * Math.SQRT2, 5);
  expect(hi).toBeCloseTo(2 * Math.SQRT2, 5);
  const before = await draw.locator("canvas").screenshot();
  // the slider rounds to its step
  const at = await slider.evaluate((el: HTMLInputElement) => {
    el.value = "0.8";
    el.dispatchEvent(new Event("input", { bubbles: true }));
    return Number(el.value).toFixed(3);
  });
  expect(Number(at)).toBeCloseTo(0.8, 1);
  await expect(page.locator('[data-nd="readout"]')).toHaveText(`x4 = ${at}`);
  await expect(draw).toHaveAttribute("aria-label", new RegExp(`cut by x4 = ${at.replace(".", "\\.")}`));
  expect((await draw.locator("canvas").screenshot()).equals(before)).toBe(false);
  // another axis: back at the centre, and the slider's range is that axis's
  await page.selectOption('[data-nd="axis"]', "1");
  await expect(page.locator('[data-nd="readout"]')).toHaveText("x2 = 0.000");
  await expect(draw).toHaveAttribute("aria-label", /cut by x2 = 0\.000/);
  // beyond the container the slice is empty and nothing is drawn, without an error
  await slider.evaluate((el: HTMLInputElement) => {
    el.value = el.max;
    el.dispatchEvent(new Event("input", { bubbles: true }));
  });
  await expect(draw).toHaveAttribute("aria-label", /cut by x2 = 2\.828/);
  await page.getByRole("button", { name: "Centre" }).click();
  await expect(page.locator('[data-nd="readout"]')).toHaveText("x2 = 0.000");

  // the drawing turns with the keyboard
  const still = await draw.locator("canvas").screenshot();
  await draw.focus();
  await page.keyboard.press("Shift+ArrowRight");
  expect((await draw.locator("canvas").screenshot()).equals(still)).toBe(false);

  // the table of the pieces is under the drawing, closed
  const pieces = page.locator(".nd-details");
  await expect(pieces).not.toHaveAttribute("open", "");
  await pieces.locator("summary").click();
  await expect(pieces.locator("tbody tr")).toHaveCount(12);

  // another P of the same shapes keeps the axis and the view
  await page.selectOption("#space-f-P", "11");
  await shown(page, "ndim/4/cube/ortho/11");
  await expect(page.locator('[data-nd="axis"]')).toHaveValue("1");
  expect(errors).toEqual([]);
});

test("5D in 5D shows the pieces in a table and a chart of two axes", async ({ page }) => {
  const errors = pageErrors(page);
  await page.goto("/");
  await page.goto("/#space=ndim/5/cube/ortho/12");
  await page.reload();
  await shown(page, "ndim/5/cube/ortho/12", LAZY_VIEW);
  const chart = page.locator(".nd-chart");
  await expect(chart.locator("g[data-piece]")).toHaveCount(12, LAZY_VIEW);
  await expect(page.locator(".nd-draw")).toHaveAttribute("aria-label", /^12 hypercubes in an orthoplex of scale \d+\.\d{6}, projected onto x1 and x2\.$/);
  const table = page.locator(".nd-table");
  await expect(table.locator("thead th")).toHaveText(["#", "Piece", "Centre", "Turn"]);
  await expect(table.locator("tbody tr")).toHaveCount(12);
  await expect(table.locator("tbody tr").first().locator("td").nth(2)).toHaveText(/^\((−?\d+\.\d{4}, ){4}−?\d+\.\d{4}\)$/);
  await expect(table.locator("tbody tr").first().locator("td").nth(3)).toHaveText(/^(not turned|axes permuted|turned about \d+°)$/);
  // other axes draw another picture, and the two axes never coincide
  const first = () => chart.locator('g[data-piece="1"] > *').first().evaluate((e) => e.outerHTML);
  const a = await first();
  await page.selectOption('[data-nd="axis-x"]', "2");
  await page.selectOption('[data-nd="axis-y"]', "4");
  await expect(page.locator(".nd-draw")).toHaveAttribute("aria-label", /projected onto x3 and x5\.$/);
  expect(await first()).not.toBe(a);
  await page.selectOption('[data-nd="axis-x"]', "4");
  await expect(page.locator('[data-nd="axis-y"]')).not.toHaveValue("4");
  // a larger dimension: 10D balls in a ball, drawn as disks
  await page.goto("/#space=ndim/10/ball/ball/5");
  await shown(page, "ndim/10/ball/ball/5");
  await runFor(page, 300);
  await expect(page.locator(".nd-chart g[data-piece] circle")).toHaveCount(5);
  await expect(page.locator(".nd-table tbody tr").first().locator("td").nth(3)).toHaveText("none (a ball)");
  expect(errors).toEqual([]);
});

test("a run on an nD problem without a record, and a submission, work through the worker and the server", async ({ page }) => {
  const errors = pageErrors(page);
  await page.goto("/");
  await page.goto("/#space=ndim/5/ball/ortho/4");
  await page.reload();
  await shown(page, "ndim/5/ball/ortho/4", LAZY_VIEW);
  // a problem without a research record: nothing to draw, a new arrangement to start from
  await expect(page.locator("#space-empty")).toBeVisible();
  await expect(page.locator("#space-solve-start")).toHaveValue("fresh");
  await expect(page.locator(".nd-chart g[data-piece]")).toHaveCount(0);
  const status = await runFor(page, 300);
  expect(status).toMatch(/beats the problem's best/);
  await expect(page.locator(".nd-chart g[data-piece]")).toHaveCount(4);
  await expect(page.locator(".nd-table tbody tr")).toHaveCount(4);
  await expect(page.locator("#space-shown")).toHaveText("the optimizer's result");
  const result = await page.evaluate(
    () => (window as unknown as { polypack: { space: { lastResult: { score: number; valid: boolean; submittable: boolean } } } }).polypack.space.lastResult,
  );
  expect(result.valid).toBe(true);
  expect(result.submittable).toBe(true);

  // the server's answer to a check it queues is shown: first its wait, then a full queue (503) and nothing stored
  await page.route("**/api/spaces/submit", async (route) => {
    await new Promise((r) => setTimeout(r, 400));
    await route.fulfill({
      status: 503,
      contentType: "application/json",
      body: JSON.stringify({ valid: false, reasons: ["the server is busy checking other entries; try again in a few seconds"], score: null }),
    });
  });
  await page.fill("#name", "ndim tester");
  await page.click("#space-submit");
  await expect(page.locator("#space-submit-result")).toHaveText("Submitting…");
  await expect(page.locator("#space-submit-result")).toHaveText(
    "✗ The server is busy checking other entries; try again in a few seconds. Nothing was stored; try again later.",
  );
  await page.unroute("**/api/spaces/submit");

  // the real server keeps it
  await page.click("#space-submit");
  await expect(page.locator("#space-submit-result")).toContainText("✓ Kept. Server score", { timeout: 30_000 });
  const api = (await (await page.request.get("/api/spaces/ndim/problems/ndim/5/ball/ortho/4")).json()) as SpaceProblemResponse;
  expect(api.site_best?.name).toBe("ndim tester");
  // the only console error is the browser's note of the 503 that the test sent itself
  expect(errors.filter((e) => !/status of 503/.test(e))).toEqual([]);
});

test("4D on 5D to 9D on 10D show the centres per facet in a chart and a table", async ({ page }) => {
  const errors = pageErrors(page);
  await page.goto("/");
  await page.goto("/#space=boundary/6cube/ball/9");
  await page.reload();
  await shown(page, "boundary/6cube/ball/9", LAZY_VIEW);
  // the d-cubes have no research record: a run gives the first packing
  await expect(page.locator("#space-empty")).toBeVisible();
  await runFor(page, 2000);
  const chart = page.locator(".nd-chart");
  await expect(chart.locator("g[data-piece]")).toHaveCount(9, LAZY_VIEW);
  await expect(page.locator(".nd-facets")).toContainText(/^Pieces per facet: [+−]x\d \d+/);
  const table = page.locator(".nd-table");
  await expect(table.locator("thead th")).toHaveText(["#", "Facet", "Centre in the facet", "Centre in x1 to x6"]);
  await expect(table.locator("tbody tr")).toHaveCount(9);
  await expect(table.locator("tbody tr").first().locator("td").nth(1)).toHaveText(/^[+−]x[1-6]$/);
  await expect(page.locator(".nd-draw")).toHaveAttribute("aria-label", /^9 balls on the boundary of the 6-cube of edge \d+\.\d{6}, their centres projected onto x\d+ and x\d+\.$/);
  // the axes are the user's
  await page.selectOption('[data-nd="axis-x"]', "3");
  await expect(page.locator(".nd-draw")).toHaveAttribute("aria-label", /projected onto x4 and x\d+\.$/);
  // cubes on the same container: a row per piece too
  await page.selectOption("#space-f-piece", "cube");
  await shown(page, "boundary/6cube/cube/9");
  await expect(page.locator("#space-empty")).toBeVisible();
  expect(errors).toEqual([]);
});

test("3D on 4D still draws all six of its containers, in the three.js views", async ({ page }) => {
  test.setTimeout(240_000);
  const errors = pageErrors(page);
  for (const polytope of ["tesseract", "5cell", "16cell", "24cell", "120cell", "600cell"]) {
    // only the tesseract has a research record; a run gives the others their first packing to draw (one ball, which
    // every container takes at its lower bound)
    const id = polytope === "tesseract" ? "boundary/tesseract/ball/3" : `boundary/${polytope}/ball/1`;
    await page.goto("/");
    await page.goto(`/#space=${id}`);
    await page.reload();
    await shown(page, id, LAZY_VIEW);
    if (polytope !== "tesseract") {
      await expect(page.locator("#space-empty")).toBeVisible();
      await runFor(page, 3000);
    }
    // the board's own 3D views, not a table: a canvas that drew something
    const canvas = page.locator(".pp3 canvas");
    await expect(canvas).toBeVisible(LAZY_VIEW);
    await expect(page.locator(".nd-chart")).toHaveCount(0);
    await expect.poll(async () => distinctColours(page, await canvas.screenshot()), { timeout: 30_000 }).toBeGreaterThan(4);
    await expect(page.locator("#family-select")).toHaveValue("boundary");
    await expect(page.locator("main")).toHaveAttribute("data-dimension", "3-on-4");
  }
  expect(errors).toEqual([]);
});

test("the Records page files the nD dimensions under their tabs and filters", async ({ page }) => {
  const errors = pageErrors(page);
  await page.goto("/records.html#space=ndim&d=6");
  const tabs = page.getByRole("navigation", { name: "Dimension" });
  await expect(tabs.locator("a")).toHaveCount(DIMENSIONS.length);
  await expect(tabs.locator('a[aria-current="page"]')).toHaveText("6D in 6D");
  await expect(page.locator("#pick-d option")).toHaveText(["6D"]);
  await expect(page.locator("#pick-piece option").first()).toHaveText("Simplex");
  await expect(page.locator("#pick-container option").last()).toHaveText("Ball");
  await expect(page.locator("#kind-summary")).toContainText("Equal 6D polytopes or balls");
  await expect(page).toHaveTitle(/Records: 6D in 6D, Polytopes and balls/);
  // a family with a record: its table, with the source of the records written under the pickers
  await page.selectOption("#pick-piece", "cube");
  await page.selectOption("#pick-container", "ortho");
  await expect(page.locator("#chart-title")).toContainText("6D: Hypercubes in an orthoplex");
  await expect(page.locator("#source-note")).toContainText("polyhedrainpolyhedra.wiki");
  await expect(page.locator("#source-note")).toContainText("Zarzuelo");
  await expect(page.locator("#series-note")).toContainText("of this family's 23 problems");
  await expect(page.locator("#rec-rows tr").first().getByRole("link", { name: /^Board/ })).toHaveAttribute("href", /^\/#space=ndim\/6\/cube\/ortho\/\d+$/);
  // all of the dimension's records, with columns that name the shapes
  await page.locator('input[name="scope"][value="all"]').check();
  await expect(page.locator("#rec-head th").first()).toHaveText("Dimension");
  await expect(page.locator("#rec-rows tr").first().locator("td").first()).toHaveText("6D");
  const ids = await page.locator("#rec-rows tr").evaluateAll((trs) => trs.map((tr) => (tr as HTMLElement).dataset.id!));
  expect(ids.length).toBeGreaterThan(100);
  expect(ids.every((id) => id.startsWith("ndim/6/"))).toBe(true);

  // 7D has no record yet, and says so
  await tabs.getByRole("link", { name: "7D in 7D" }).click();
  await expect(tabs.locator('a[aria-current="page"]')).toHaveText("7D in 7D");
  await expect(page.locator("#source-note")).toHaveText(/^No research record yet/);
  await expect(page).toHaveURL(/#space=ndim&d=7/);

  // the d-cube tab, and the older address of 3D on 4D
  await tabs.getByRole("link", { name: "4D on 5D" }).click();
  await expect(tabs.locator('a[aria-current="page"]')).toHaveText("4D on 5D");
  await expect(page.locator("#pick-polytope option")).toHaveText(["5-cube"]);
  await page.goto("/records.html#space=boundary&polytope=tesseract");
  await expect(tabs.locator('a[aria-current="page"]')).toHaveText("3D on 4D");
  await expect(page.locator("#pick-polytope option")).toHaveCount(6);
  await page.goto("/records.html#space=boundary");
  await expect(tabs.locator('a[aria-current="page"]')).toHaveText("3D on 4D");
  expect(errors).toEqual([]);
});

test("the dimensions page has a card for each nD dimension and links it to the board and the records", async ({ page }) => {
  await page.goto("/dimensions.html");
  const sections = page.locator("#dimension-groups > section");
  await expect(sections).toHaveCount(DIMENSIONS.length);
  const four = page.locator('section[data-dimension="4-in-4"]');
  await expect(four.locator("h3")).toHaveText("4D in 4D");
  await expect(four.locator("p.dimension-summary")).toHaveText(/^Equal 4D polytopes or balls/);
  await expect(four.locator("li h4")).toHaveText("Polytopes and balls");
  await expect(four.getByRole("link", { name: /^Board/ })).toHaveAttribute("href", "/#space=ndim/4/cube/cube/12");
  await expect(four.getByRole("link", { name: /^Records/ })).toHaveAttribute("href", "#space=ndim&d=4");
  await expect(four).toContainText("1,911");
  const cube = page.locator('section[data-dimension="9-on-10"]');
  await expect(cube.getByRole("link", { name: /^Board/ })).toHaveAttribute("href", "/#space=boundary/10cube/ball/9");
  await expect(cube).toContainText("800");
  // a keyed card opens its own dimension on the board
  await four.getByRole("link", { name: /^Board/ }).click();
  await shown(page, "ndim/4/cube/cube/12", LAZY_VIEW);
  expect(DIMENSIONS.filter((d) => d.families.some((f) => familyKey(f).startsWith("ndim-"))).length).toBe(7);
});

for (const scheme of ["light", "dark"] as const) {
  test(`the nD drawings and their tables have no automatic accessibility violations (${scheme})`, async ({ browser }) => {
    const page = await (await browser.newContext({ colorScheme: scheme })).newPage();
    for (const [id, ready] of [
      ["ndim/4/cube/ortho/12", ".nd-draw canvas"],
      ["ndim/5/cube/ortho/12", ".nd-chart g[data-piece]"],
    ] as const) {
      await page.goto(`/#space=${id}`);
      await expect(page.locator(ready).first()).toBeVisible(LAZY_VIEW);
      if (id.startsWith("ndim/4")) await page.locator(".nd-details summary").click();
      await expect(page.locator(".nd-table tbody tr").first()).toBeVisible();
      const results = await new AxeBuilder({ page }).analyze();
      expect(results.violations.map((v) => `${v.id}: ${v.nodes.map((n) => n.target.join(" ")).join(", ")}`), id).toEqual([]);
    }
    await page.close();
  });
}
