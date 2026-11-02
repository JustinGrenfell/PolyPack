import { expect, test } from "@playwright/test";
import { LAZY_VIEW } from "./helpers";

// The nD dimensions at phone width: the page does not scroll sideways, the drawing is as tall as the screen is wide and its
// table scrolls in its own box.

const overflow = (page: import("@playwright/test").Page) =>
  page.evaluate(() => document.documentElement.scrollWidth - document.documentElement.clientWidth);

test("phone: 4D in 4D, 5D in 5D and a d-cube fit, with their drawing and table", async ({ page }) => {
  await page.goto("/#space=ndim/4/cube/ortho/12");
  await expect(page.locator("#space-problem-note")).toHaveText("Problem ndim/4/cube/ortho/12", LAZY_VIEW);
  const canvas = page.locator(".nd-draw canvas");
  await expect(canvas).toBeVisible(LAZY_VIEW);
  const width = page.viewportSize()!.width;
  const box = (await canvas.boundingBox())!;
  expect(box.width).toBeLessThanOrEqual(width + 1);
  // as tall as the screen is wide, at most 70% of its height
  expect(box.height).toBeGreaterThan(width * 0.6);
  expect(await overflow(page)).toBeLessThanOrEqual(1);
  // the slider and the axis are reachable by touch
  const slider = page.locator('[data-nd="position"]');
  await slider.scrollIntoViewIfNeeded();
  const s = (await slider.boundingBox())!;
  expect(s.x).toBeGreaterThanOrEqual(0);
  expect(s.x + s.width).toBeLessThanOrEqual(width + 1);
  expect(s.width).toBeGreaterThan(100);
  await page.selectOption('[data-nd="axis"]', "0");
  await expect(page.locator(".nd-draw")).toHaveAttribute("aria-label", /cut by x1 = 0\.000/);

  await page.goto("about:blank");
  await page.goto("/#space=ndim/8/ortho/cube/6");
  await expect(page.locator("#space-problem-note")).toHaveText("Problem ndim/8/ortho/cube/6", LAZY_VIEW);
  await page.locator("#space-empty").waitFor({ state: "attached" });
  expect(await overflow(page)).toBeLessThanOrEqual(1);

  // a 5D record: the chart is as wide as the screen, the table's box scrolls on its own
  await page.goto("about:blank");
  await page.goto("/#space=ndim/5/cube/ortho/12");
  await expect(page.locator("#space-problem-note")).toHaveText("Problem ndim/5/cube/ortho/12", LAZY_VIEW);
  await expect(page.locator(".nd-chart g[data-piece]")).toHaveCount(12, LAZY_VIEW);
  const chart = (await page.locator(".nd-chart").boundingBox())!;
  expect(chart.width).toBeLessThanOrEqual(width + 1);
  expect(chart.width).toBeGreaterThan(width * 0.8);
  const wrap = page.locator(".nd-table-wrap");
  await wrap.scrollIntoViewIfNeeded();
  expect(await wrap.evaluate((e) => e.scrollWidth > e.clientWidth)).toBe(true);
  expect(await overflow(page)).toBeLessThanOrEqual(1);
  await expect(page.locator(".nd-table tbody tr")).toHaveCount(12);

  // 10D: the centre has ten numbers, so the table scrolls further and the page still does not
  await page.goto("about:blank");
  await page.goto("/#space=ndim/10/ball/ball/5");
  await expect(page.locator("#space-problem-note")).toHaveText("Problem ndim/10/ball/ball/5", LAZY_VIEW);
  expect(await overflow(page)).toBeLessThanOrEqual(1);

  // the menu names the dimension in full
  const menu = page.locator("#family-select");
  await menu.scrollIntoViewIfNeeded();
  await expect(menu.locator("option:checked")).toHaveText("10D in 10D: Polytopes and balls");
  const m = (await menu.boundingBox())!;
  expect(m.x + m.width).toBeLessThanOrEqual(width + 1);
});

test("phone: the records of an nD dimension fit", async ({ page }) => {
  await page.goto("/records.html#space=ndim&d=5");
  await expect(page.locator('.dim-tabs a[aria-current="page"]')).toHaveText("5D in 5D");
  await expect(page.locator("#series-note")).toContainText("this family's 31 problems");
  expect(await overflow(page)).toBeLessThanOrEqual(1);
  await page.goto("/records.html#space=boundary&polytope=6cube");
  await expect(page.locator('.dim-tabs a[aria-current="page"]')).toHaveText("5D on 6D");
  expect(await overflow(page)).toBeLessThanOrEqual(1);
});
