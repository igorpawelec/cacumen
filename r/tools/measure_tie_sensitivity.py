"""How much do the crowns move when only the plateau ties change?

A canopy height model is full of plateaus — flat roofs of smoothed crowns,
and whole regions of ground clipped to the same value. Tree detection and
marker-based watershed both have to break those ties somehow, and the choice
is arbitrary: scikit-image floods first-in-first-out by insertion age, a
different implementation may not, and neither is more correct than the other.

This script puts a number on what that arbitrariness costs. It adds noise
far below any physical or instrumental meaning, which changes nothing about
the canopy but does break the ties differently each time. Whatever the crowns
do in response is caused by tie-breaking alone. Two levels of noise, because
the stages of the method work at two precisions:

(a) 1e-9 m, a nanometre, with the tops free to move. Tree detection runs in
    double precision and sees it. The region grower stores its surface in
    single precision, where a nanometre on a height of metres is nothing: it
    is handed the same surface every time. What moves here moves because
    the tops did.

(b) 1e-5 m, ten micrometres, with the tops held where the undisturbed run
    found them. This is a few units in the last place of a single-precision
    height, so it reaches the watershed and the region growing; with the
    tops frozen, what moves here moves because of them alone.

Until 2026-11 branch (b) was run with the nanometre as well, and reported
that nothing moved. Nothing could: the grower never saw the noise. The count
of pixels whose single-precision height differs from the undisturbed run is
printed for both levels, so that a measurement of nothing shows as one.

    python3 tools/measure_tie_sensitivity.py [n_repeats]

Copyright (C) 2025 Igor Pawelec. Licence: GPLv3.
"""
import os
import sys

import numpy as np
import rasterio
from pycacumen import (smooth_chm, detect_tops, as_pixels,
                   HierarchicalRegionGrower)

SRC = os.environ.get("CHM_DIR", os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                             "..", "..", "python", "test_data"))
N = int(sys.argv[1]) if len(sys.argv) > 1 else 8
EPS = 1e-9          # (a): seen by detection only
EPS_GROW = 1e-5     # (b): survives the grower's single precision

HMIN, WS, SMOOTH_WS = 5.0, 5, 3
VT, MASK, RULE = 8.0, 1.0, "height"


def grow(sm, tops):
    g = HierarchicalRegionGrower(sm)
    return g.run_all(tops, variance_thresh=VT, mask_thresh=MASK,
                     conflict_rule=RULE)


def run(sm):
    tops = as_pixels(detect_tops(sm, hmin=HMIN, ws=WS))
    if not tops:
        return None, None, 0
    return grow(sm, tops), tops, len(tops)


def seed_map(out, tops):
    """Crown ids -> the tree top each pixel belongs to.

    Crown ids are arbitrary: id 7 in one run and id 7 in another need not
    be the same tree. Comparing the integers directly counts a renumbering
    as a total change and reports ~100% wherever the top ordering shifts,
    which is a measurement artefact, not a result. Resolving each pixel to
    its seed *coordinates* compares the partition itself.
    """
    coords = np.zeros((len(tops) + 1, 2), dtype=np.int64)
    coords[0] = (-1, -1)                       # background
    for i, (r, c) in enumerate(tops, start=1):
        coords[i] = (r, c)
    idx = np.clip(out, 0, len(tops))
    return coords[idx]


def seen_by_grower(sm, base_sm):
    """Masked pixels whose height differs from the baseline in single precision."""
    return int(((sm.astype(np.float32) != base_sm.astype(np.float32))
                & (base_sm > MASK)).sum())


print("%d repeats per scene; (a) noise %g m, tops free; (b) noise %g m, tops fixed\n"
      % (N, EPS, EPS_GROW))
print("%-16s %6s %10s %9s %12s %17s %17s" %
      ("case", "tops", "canopy px", "tops mvd", "tops +/-",
       "(a) worst", "(b) worst"))
print("-" * 97)

rows = []
seen_a = seen_b = 0
for f in sorted(os.listdir(SRC)):
    if not f.endswith(".tif"):
        continue
    name = f[:-4]
    with rasterio.open(os.path.join(SRC, f)) as src:
        raw = src.read(1).astype(np.float64)
        nod = src.nodata
    if nod is not None:
        raw[raw == nod] = 0.0
    raw[~np.isfinite(raw)] = 0.0

    base_sm = smooth_chm(raw, ws=SMOOTH_WS, method="median")
    base, base_tops, n_tops = run(base_sm)
    if base is None:
        continue
    canopy = int((base > 0).sum())
    base_seed = seed_map(base, base_tops)

    diffs, fixed_diffs, top_moves, top_counts = [], [], [], []
    for k in range(N):
        rng = np.random.default_rng(1000 + k)
        unit = rng.uniform(0, 1, raw.shape)

        # (a) a nanometre, whole pipeline: the tops are free to move, because
        #     a local maximum on a plateau is itself a tie.
        sm = smooth_chm(raw + EPS * unit, ws=SMOOTH_WS, method="median")
        seen_a += seen_by_grower(sm, base_sm)
        out, tops, n_k = run(sm)
        if out is None:
            continue
        # How many tops are not where the undisturbed run has them (gone,
        # new or moved, each counted once per position), and how the number
        # of tops itself changes.
        top_moves.append(len(set(tops) ^ set(base_tops)))
        top_counts.append(n_k - n_tops)
        moved = (seed_map(out, tops) != base_seed).any(axis=-1)
        # Only canopy counts: background staying background is not news.
        diffs.append(int((moved & ((out > 0) | (base > 0))).sum()))

        # (b) ten micrometres, tops frozen at the baseline, so detection is
        #     out of it and what is left is the watershed and the region
        #     growing alone. Without this split the two effects are reported
        #     as one and the watershed gets blamed for what detection did.
        sm_g = smooth_chm(raw + EPS_GROW * unit, ws=SMOOTH_WS, method="median")
        seen_b += seen_by_grower(sm_g, base_sm)
        out_f = grow(sm_g, base_tops)
        moved_f = (seed_map(out_f, base_tops) != base_seed).any(axis=-1)
        fixed_diffs.append(int((moved_f & ((out_f > 0) | (base > 0))).sum()))

    if not diffs:
        continue
    d = np.array(diffs, float)
    fd = np.array(fixed_diffs, float) if fixed_diffs else np.zeros(1)
    print("%-16s %6d %10d %9.1f %+5d..%+-5d %9d %6.2f%% %9d %6.2f%%" %
          (name, n_tops, canopy, float(np.mean(top_moves)),
           min(top_counts), max(top_counts),
           int(d.max()), 100 * d.max() / max(canopy, 1),
           int(fd.max()), 100 * fd.max() / max(canopy, 1)))
    rows.append((name, n_tops, canopy, d.max(), fd.max()))

if rows:
    tot_canopy = sum(r[2] for r in rows)
    tot_free = sum(r[3] for r in rows)
    tot_fixed = sum(r[4] for r in rows)
    print("-" * 97)
    print("across %d scenes, worst case of %s canopy pixels:" %
          (len(rows), format(tot_canopy, ",")))
    print("  (a) %g m, tops free  : %d change crown (%.2f%%)" %
          (EPS, int(tot_free), 100 * tot_free / tot_canopy))
    print("  (b) %g m, tops fixed : %d change crown (%.2f%%)" %
          (EPS_GROW, int(tot_fixed), 100 * tot_fixed / tot_canopy))
    print("masked pixels the grower saw change, all repeats: (a) %s, (b) %s" %
          (format(seen_a, ","), format(seen_b, ",")))
    print("\nNothing about the canopy changed. Only which pixel won a tie.")
