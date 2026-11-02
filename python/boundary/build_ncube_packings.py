"""Test packings on the boundary of the d-cube, d = 5..10 (no records: BOUNDARY.md 4B), certified exactly by ncube.py.

Writes fixtures/ncube_packings.json: for each d a few small constructions (the engine and the server tests read them too), each
scaled by lam = 1 + 1e-9 (the shrink argument of BOUNDARY.md 4.4: a valid layout at s becomes strictly valid at lam s with its
centres scaled by lam) and certified strictly valid. The big tilings (P = 2d k^(d-1) cubes at s = k, k^n balls per cell) are built by the
tests themselves. Run: python3 python/boundary/build_ncube_packings.py
"""

import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import certify as C  # noqa: E402
import ncube as N  # noqa: E402

LAM = 1 + 1e-9


def scaled(s0, pieces):
    return {"s": s0 * LAM, "pieces": [{"cell": c, "u": [x * LAM for x in u]} for c, u in pieces]}


def ball_centre1(d):
    return scaled(math.sqrt(0.5), [(0, [0.0] * (d - 1))])


def ball_cells(d):
    return scaled(1.0, [(c, [0.0] * (d - 1)) for c in range(2 * d)])


def ball_diag2(d):
    n = d - 1
    s0 = 1 + 1 / math.sqrt(n)
    a = (s0 - 1) / 2
    pieces = []
    for c in range(2 * d):
        p = [(-a if (c >> k) & 1 else a) for k in range(n)]
        pieces += [(c, p), (c, [-x for x in p])]
    return scaled(s0, pieces)


def cube_centre1(d):
    return scaled(1.0, [(0, [0.0] * (d - 1))])


def cube_cells(d):
    return scaled(1.0, [(c, [0.0] * (d - 1)) for c in range(2 * d)])


def ring_positions(d, s, count, phase, ra=0, rb=1):
    """`count` points of the ring of the cells (+-e_ra), (+-e_rb): [(cell, u)] with the other coordinates 0 (cross section centre)."""
    h = s / 2
    out = []
    for j in range(count):
        t = (phase + j * 4 * s / count) % (4 * s)
        tt = (t + h) % (4 * s)
        seg = int(tt // s)
        off = tt - seg * s - h
        x = [0.0] * d
        if seg == 0:
            x[ra], x[rb], c = h, off, 2 * ra + 1
        elif seg == 1:
            x[rb], x[ra], c = h, -off, 2 * rb + 1
        elif seg == 2:
            x[ra], x[rb], c = -h, -off, 2 * ra
        else:
            x[rb], x[ra], c = -h, off, 2 * rb
        ax = N.axes(d, c >> 1)
        out.append((c, [x[a] for a in ax]))
    return out


def cube_ring_ridges(d):
    """4 unit cubes centred on the 4 ridges of the ring (+-e_0, +-e_1) at s = 1: each crosses one axis and lies half in two cells."""
    return scaled(1.0, ring_positions(d, 1.0, 4, 0.5))


def cube_ring_shifted(d):
    """8 cubes round the ring at s = 2 (spacing 1), half a step off the cells' centres: some straddle ridges."""
    return scaled(2.0, ring_positions(d, 2.0, 8, 0.5))


CASES = [
    ("ball", "centre1", ball_centre1, lambda d: 1),
    ("ball", "cells", ball_cells, lambda d: 2 * d),
    ("ball", "diag2", ball_diag2, lambda d: 4 * d),
    ("cube", "centre1", cube_centre1, lambda d: 1),
    ("cube", "cells", cube_cells, lambda d: 2 * d),
    ("cube", "ring_ridges", cube_ring_ridges, lambda d: 4),
    ("cube", "ring_shifted", cube_ring_shifted, lambda d: 8),
]


def main():
    packings = []
    for d in range(5, 11):
        for piece, kind, make, count in CASES:
            cfg = make(d)
            P = count(d)
            assert len(cfg["pieces"]) == P
            pid = f"boundary/{d}cube/{piece}/{P}"
            cert = C.certify(pid, cfg)
            packings.append({"id": pid, "polytope": f"{d}cube", "piece": piece, "P": P, "kind": kind, "score_cert": cfg["s"], "config": cfg, "margins": {k: v for k, v in cert.items() if k.startswith("min_")}})
            print(pid, kind, cfg["s"], cert)
    with open(os.path.join(HERE, "fixtures", "ncube_packings.json"), "w") as fh:
        json.dump({"note": "built by python/boundary/build_ncube_packings.py; every packing certified strictly valid by ncube.py", "packings": packings}, fh)


if __name__ == "__main__":
    main()
