"""Builds python/boundary/fixtures/platonic_packings.json: the test packings of Platonic solids on the regular 4-polytopes (no records;
BOUNDARY.md 4C.6), each certified strictly valid by the exact certifier platonic.py.

Input: fixtures/platonic_packings_raw.json, written by packages/boundary/scripts/gen-platonic-packings.ts (the cells themselves at s = 1,
and short runs of the boundary search: touching, valid at tolerance 0). Each is scaled about the centres of the cells as BOUNDARY.md section 4.4
does for the records: lam_k = 1 + k 1e-12 for k = 0, 1, 2, ..., the first at which the scaled binary64 configuration certifies strictly valid, also with twice the inflation.
Strict validity is decided on the solids inflated by the factor 1 + EPS_STRICT (platonic.py), which exceeds the reading discrepancy
many times over; the stored margin is the largest power of two eps, up to 2^-30, for which the inflated solids are still valid.

    python3 python/boundary/build_platonic_packings.py
"""
import json
import os
import sys
import time
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import certify as C  # noqa: E402
import platonic  # noqa: E402

RAW = os.path.join(HERE, "fixtures", "platonic_packings_raw.json")
OUT = os.path.join(HERE, "fixtures", "platonic_packings.json")


def scaled(config, lam):
    out = {"s": config["s"] * lam, "pieces": []}
    for p in config["pieces"]:
        q = {"cell": p["cell"], "u": [x * lam for x in p["u"]]}
        if "q" in p:
            q["q"] = p["q"]
        out["pieces"].append(q)
    return out


def margin_eps(polytope, piece, P, cfg):
    """The largest eps = 2^-j (j = 39 down to 30) at which the inflated solids are still valid; EPS_STRICT = 2^-40 itself otherwise."""
    best = platonic.EPS_STRICT
    for j in range(39, 29, -1):
        eps = Fraction(1, 2 ** j)
        try:
            v = platonic.check_pieces(polytope, piece, P, cfg, C.to_fraction, strict=True, eps=eps)
        except platonic.PlatonicError:
            break
        if not v.ok:
            break
        best = eps
    return best


def build(raw, log=print):
    out = []
    for r in raw["packings"]:
        t = time.time()
        polytope, piece, P = r["polytope"], r["piece"], r["P"]
        for k in range(0, 400):
            lam = 1 + k * 1e-12
            cfg = scaled(r["config"], lam) if k else r["config"]
            v = C.check(r["id"], cfg, strict=True)
            # a margin of at least twice the inflation the certifier asks (the engine's reading is far inside the first)
            if v.ok and platonic.check_pieces(polytope, piece, P, cfg, C.to_fraction, strict=True, eps=2 * platonic.EPS_STRICT).ok:
                break
        else:
            raise SystemExit(f"{r['id']}: no strict scaling found")
        summ = C.summary(v)
        eps = margin_eps(polytope, piece, P, cfg)
        disc = platonic.reading_discrepancy(polytope, piece, cfg, C.to_fraction)
        out.append(
            {
                "id": r["id"],
                "polytope": polytope,
                "piece": piece,
                "P": P,
                "source": r["source"],
                "score_cert": cfg["s"],
                "k": k,
                "config": cfg,
                "inflation": float(platonic.EPS_STRICT),
                "margin_inflation": float(eps),
                "reading_discrepancy": disc,
                "chains": summ["chains"],
                "parts": summ["parts"],
            }
        )
        log(f"{r['id']} ({r['source']}): s = {cfg['s']!r} (k = {k}), margin eps {float(eps):.3g}, reading {disc:.2g}, {time.time() - t:.1f} s")
    return out


def main():
    with open(RAW) as fh:
        raw = json.load(fh)
    packings = build(raw)
    with open(OUT, "w") as fh:
        json.dump({"note": "built by python/boundary/build_platonic_packings.py from platonic_packings_raw.json; every packing certified strictly valid by platonic.py", "packings": packings}, fh)
    print(f"wrote {len(packings)} packings")


if __name__ == "__main__":
    main()
