"""Builds python/boundary/fixtures/nd_packings.json: the test packings of balls on the simplex, orthoplex and demicube, d = 5..10
(no records; BOUNDARY.md 4D.7), each certified strictly valid by the exact certifier ndpoly.py.

Input: fixtures/nd_packings_raw.json, written by packages/boundary/scripts/gen-nd-packings.ts (exact layouts and short seeded
searches, drawn in by the kernel's finish(): touching, valid at tolerance 0). Each is scaled about the centres of the cells as BOUNDARY.md
section 4.4 does for the records: lam_k = 1 + k 1e-12 for k = 1, 2, ..., the first at which the scaled binary64 configuration
certifies strictly valid with every margin (clearance over 1/2, pair distance and loop length over 1) at least MIN_MARGIN = 1e-12:
the engine's floats are read to within about 1e-16 s of the point they mean (BOUNDARY.md 4D.4), and the strict test keeps a slack for that,
so smaller margins would not bind the configuration as the engine reads it. The result is stored with its exact margins (rounded down).

    python3 python/boundary/build_nd_packings.py
"""
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import certify as C  # noqa: E402
import ndpoly  # noqa: E402
from fractions import Fraction  # noqa: E402

RAW = os.path.join(HERE, "fixtures", "nd_packings_raw.json")
OUT = os.path.join(HERE, "fixtures", "nd_packings.json")


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
            # extra = 2e-12: the search looks for margins up to that, so that a margin below MIN_MARGIN is seen and one that is not seen is above it
            v = ndpoly.check(r["polytope"], "ball", r["P"], cfg, C.to_fraction, strict=True, extra=Fraction(2, 10 ** 12))
            if v.ok and margins_ok(ndpoly.summary(v)):
                break
        else:
            raise SystemExit(f"{r['id']}: no strict scaling found")
        summ = ndpoly.summary(v)
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
        json.dump({"note": "built by python/boundary/build_nd_packings.py from nd_packings_raw.json; every packing certified strictly valid by ndpoly.py", "packings": packings}, fh)
    print(f"wrote {len(packings)} packings")


if __name__ == "__main__":
    main()
