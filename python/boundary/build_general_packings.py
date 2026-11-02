"""Builds python/boundary/fixtures/general_packings.json: the test packings of balls on the five containers beyond the tesseract
(no records; BOUNDARY.md 4A.7), each certified strictly valid by the exact certifier.

Input: fixtures/general_packings_raw.json, written by packages/boundary/scripts/gen-general-packings.ts (random sequential layout
drawn in by the engine's finish(): touching, valid at tolerance 0). Each is scaled about the centres of the cells as BOUNDARY.md
section 4.4 does for the records: lam_k = 1 + k 1e-12 for k = 1, 2, ..., the first at which the scaled binary64 configuration
certifies strictly valid with every margin (clearance over 1/2, pair distance and loop length over 1) at least MIN_MARGIN = 1e-12:
the engine's floats are read to within 3e-15 s of the point they mean (BOUNDARY.md 4A.3), so smaller margins would not bind the
configuration as the engine reads it. The result is stored with its exact margins (rounded down).

    python3 python/boundary/build_general_packings.py
"""
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import certify as C  # noqa: E402

RAW = os.path.join(HERE, "fixtures", "general_packings_raw.json")
OUT = os.path.join(HERE, "fixtures", "general_packings.json")


MIN_MARGIN = 1e-12


def margins_ok(summ):
    """Every margin of the summary (clearance - 1/2, distance - 1, loop - 1; None when there is none) is at least MIN_MARGIN."""
    for key, limit in (("min_clearance", 0.5), ("min_distance", 1.0), ("min_loop", 1.0)):
        x = summ[key]
        if x is not None and x - limit < MIN_MARGIN:
            return False
    return True


def scaled(config, lam):
    return {"s": config["s"] * lam, "pieces": [{"cell": p["cell"], "u": [x * lam for x in p["u"]]} for p in config["pieces"]]}


def build(raw, log=print):
    out = []
    for r in raw["packings"]:
        t = time.time()
        for k in range(0, 2000):
            lam = 1 + k * 1e-12
            cfg = scaled(r["config"], lam) if k else r["config"]
            v = C.check(r["id"], cfg, strict=True)
            if v.ok and margins_ok(C.summary(v)):
                break
        else:
            raise SystemExit(f"{r['id']}: no strict scaling found")
        summ = C.summary(v)
        out.append(
            {
                "id": r["id"],
                "polytope": r["polytope"],
                "P": r["P"],
                "source": r["source"],
                "score_cert": cfg["s"],
                "k": k,
                "config": cfg,
                "min_clearance": summ["min_clearance"],
                "min_distance": summ["min_distance"],
                "min_loop": summ["min_loop"],
            }
        )
        log(f"{r['id']}: s = {cfg['s']!r} (k = {k}) in {time.time() - t:.1f} s")
    return out


def main():
    with open(RAW) as fh:
        raw = json.load(fh)
    packings = build(raw)
    with open(OUT, "w") as fh:
        json.dump({"note": "built by python/boundary/build_general_packings.py from general_packings_raw.json; every packing certified strictly valid by general.py", "packings": packings}, fh)
    print(f"wrote {len(packings)} packings")


if __name__ == "__main__":
    main()
