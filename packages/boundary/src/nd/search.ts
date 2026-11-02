/**
 * The boundary search for balls on the simplex, orthoplex and demicube (d = 5..10): the structured starts of the kernel (seeds.ts) finished by its
 * exact finishing (kernel.ts finish), best kept. It is not an optimizer of records: these problem sets ship without records
 * (BOUNDARY.md 4D), and the generic methods (packages/methods) run on the kernel for anything more. It exists so that a
 * boundary job on one of these problems has an engine with the contract of searchBoundary: seeded, budgeted in evaluations
 * (a finish costs ND_FINISH_COST), stoppable, reproducible from (problem, start, seed, evals) when no clock limit is set.
 */
import type { BoundaryConfig } from "../types";
import type { BoundarySearchOptions, BoundarySearchResult } from "../search/search";
import { mulberry32 } from "../search/random";
import { ND_FINISH_COST, ndKernel } from "./kernel";
import type { NDProblem } from "./validate";

export function searchND(problem: NDProblem, opts: BoundarySearchOptions): BoundarySearchResult {
  const k = ndKernel(problem);
  const rand = mulberry32(opts.seed);
  const clock = opts.now ?? (() => performance.now() / 1000);
  const t0 = clock();
  const maxEvals = Number.isFinite(opts.budget?.evals) ? Math.max(0, Math.floor(opts.budget.evals)) : 0;
  const maxSeconds = opts.budget?.time_s;
  let evals = 0;
  let best: { config: BoundaryConfig; score: number } | null = null;
  let startScore: number | null = null;
  const trace: BoundarySearchResult["trace"] = [];
  let ended: BoundarySearchResult["ended"] = "evals";
  const take = (
    config: BoundaryConfig,
    score: number,
    source: BoundarySearchResult["trace"][number]["source"],
    stage: "given" | "coarse",
  ) => {
    if (best && !(score < best.score)) return;
    best = { config, score };
    trace.push({ evals, seconds: clock() - t0, score, source, stage });
    opts.onImprove?.({ config, score });
  };
  if (opts.start) {
    const v = k.validate(opts.start);
    if (v.valid && v.score !== null) {
      startScore = v.score;
      take(opts.start, v.score, "start", "given");
    }
  }
  const stop = () => {
    if (opts.shouldStop?.()) return "stopped" as const;
    if (maxSeconds !== undefined && clock() - t0 > maxSeconds) return "time" as const;
    if (evals + ND_FINISH_COST > maxEvals) return "evals" as const;
    return null;
  };
  for (;;) {
    const why = stop();
    if (why) {
      ended = why;
      break;
    }
    const seed = k.seeds!(rand, 1)[0];
    evals += ND_FINISH_COST;
    const f = k.finish(seed);
    if (f) take(f.config as BoundaryConfig, f.score, "random", "coarse");
    // the bound is reached: nothing smaller exists
    if (best && k.bound !== null && (best as { score: number }).score <= k.bound * (1 + 1e-12)) {
      ended = "bound";
      break;
    }
  }
  const b = best as { config: BoundaryConfig; score: number } | null;
  return {
    config: b ? b.config : null,
    score: b ? b.score : null,
    evals,
    ended,
    trace,
    startScore,
    starts: trace.length,
    stats: { coarse: evals, polish: 0, final: 0, validations: 0 },
  };
}
