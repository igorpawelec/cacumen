"""Write reference crowns from pycacumen for tools/cross_validate_against_pycacumen.R,
tools/quick_start_on_test_tiles.R and tools/real_tiles_stage_by_stage.R.

    pip install "../python[io]"      # pycacumen from this checkout, with rasterio
    python3 tools/generate_pycacumen_reference.py      # from r/

Copyright (C) 2025 Igor Pawelec. Licence: GPLv3.
"""
import glob
import os
import numpy as np
from pycacumen import smooth_chm, detect_tops, as_pixels, HierarchicalRegionGrower

OUT = os.path.join("tools", "reference")
os.makedirs(OUT, exist_ok=True)
# The files for the real tiles are written only when rasterio is installed.
# Left over from an earlier run, they would be held against rcacumen as if
# they came from the pycacumen installed now.
for stale in glob.glob(os.path.join(OUT, "real*.csv")) + glob.glob(os.path.join(OUT, "tiles.csv")):
    os.remove(stale)

rng = np.random.default_rng(11)
cases = []
# Cases 60-69 are large and crowded on purpose. The original 60 are 30-60 px
# with 2-7 trees, which is nothing like a real canopy height model: the
# divergence found on real CHMs in 2026-07 lives in plateau structure that
# only appears once crowns are dense enough to touch.
for k in range(70):
    big = k >= 60
    n = int(rng.integers(140, 200)) if big else int(rng.integers(30, 60))
    yy, xx = np.mgrid[0:n, 0:n]
    chm = np.zeros((n, n))
    n_trees = int(rng.integers(25, 60)) if big else int(rng.integers(2, 8))
    for _ in range(n_trees):
        r, c = rng.integers(4, n - 4, 2)
        chm = np.maximum(chm, rng.uniform(10, 28) *
                         np.exp(-((yy - r) ** 2 + (xx - c) ** 2)
                                / (2 * rng.uniform(3, 7) ** 2)))
    if k % 3 == 0:
        chm = np.round(chm * 2) / 2          # plateaus, where ties are common
    if k % 7 == 0:
        # Nodata holes. Until pycacumen 0.5.0 / rcacumen 0.4.0 these two disagreed
        # completely here: scipy's comparison filters are undefined on NaN,
        # so Python smoothed the hole edges to arbitrary window elements,
        # while R propagated NA and then stopped with "missing value where
        # TRUE/FALSE needed". Both now skip NaN, and this is what checks it.
        chm[2:6, 2:6] = np.nan
        chm[n - 7:n - 3, n - 9:n - 5] = np.nan
    chm = chm.astype(np.float32)

    sm = smooth_chm(chm, ws=3, method="median")
    sub = detect_tops(sm, hmin=5, ws=5)
    tops = as_pixels(sub)
    if not tops:
        continue
    vt = float(rng.choice([2.0, 8.0, 25.0]))
    rule = str(rng.choice(["height", "distance", "similarity"]))
    # protect_seeds and retry_rejected went uncovered until 2026-07-22:
    # every case ran with both off, so a divergence in either would have
    # gone unseen. They match, and this is what keeps them matching.
    protect = bool(k % 5 == 0)
    retry = bool(k % 4 == 0)

    g = HierarchicalRegionGrower(sm)
    out = g.run_all(tops, variance_thresh=vt, mask_thresh=1.0,
                    conflict_rule=rule, protect_seeds=protect,
                    retry_rejected=retry)
    np.savetxt(f"{OUT}/sm_{k}.csv", sm, delimiter=",")
    np.savetxt(f"{OUT}/tops_{k}.csv", np.array(tops) + 1, fmt="%d", delimiter=",")
    # Subpixel tops, before as_pixels. Without these the R side of
    # detect_tops is never compared against anything: the validator used to
    # read tops that Python had already detected *and* floored, so R's own
    # detect_tops -> grow_crowns handoff went unexercised.
    np.savetxt(f"{OUT}/sub_{k}.csv", np.asarray(sub, float) + 1.0,
               delimiter=",", fmt="%.17g")
    np.savetxt(f"{OUT}/out_{k}.csv", out, fmt="%d", delimiter=",")
    cases.append((k, vt, rule, g.n_contested, int(protect), int(retry)))

with open(f"{OUT}/meta.csv", "w") as f:
    for k, vt, rule, nc, ps, rr in cases:
        f.write(f"{k},{vt},{rule},{nc},{ps},{rr}\n")
print(f"{len(cases)} reference cases written to {OUT}/")

# The README quick start on the real tiles of the repository: how many tops
# and crowns pycacumen finds, for tools/quick_start_on_test_tiles.R to hold
# rcacumen against. The synthetic cases above never call merge_tops or
# screen_tops, which is where rcacumen 0.5.0 kept a top that pycacumen
# dropped and then stopped on it. Needs rasterio, which the array API above
# does not, so without it this part is skipped rather than failed.
try:
    from pycacumen import CrownDelineator
    import rasterio  # noqa: F401  (from_file needs it; fail here, not per tile)
except ImportError as e:
    print(f"real tiles skipped ({e}); install pycacumen[io] to write tiles.csv")
else:
    tiles = {}
    for folder in ("test_data", os.path.join("..", "python", "test_data")):
        for path in sorted(glob.glob(os.path.join(folder, "*.tif"))):
            tiles.setdefault(os.path.basename(path), path)
    with open(f"{OUT}/tiles.csv", "w") as f:
        f.write("tile,tops,crowns\n")
        for name, path in sorted(tiles.items()):
            cd = CrownDelineator.from_file(path, quiet=True)
            cd.smooth(ws=3).detect(hmin=7, ws=5).merge(5.0).screen(10.0)
            crowns = cd.delineate(variance_thresh=2.0)
            f.write(f"{name},{len(cd.tops)},{len(set(np.unique(crowns)) - {0})}\n")
    print(f"{len(tiles)} real tiles written to {OUT}/tiles.csv")

    # The same tiles once more, stage by stage, for
    # tools/real_tiles_stage_by_stage.R: the heights as read, and what
    # pycacumen makes of them at every step. Each step of rcacumen can then be
    # fed pycacumen's own output of the step before, which pins a difference
    # to the stage that makes it. Read and run as
    # tools/measure_tie_sensitivity.py reads and runs them, so the two tools
    # speak of the same crowns.
    HMIN, WS, SMOOTH_WS, VT, MASK, RULE = 5.0, 5, 3, 8.0, 1.0, "height"
    with open(f"{OUT}/real.csv", "w") as f:
        f.write("tile,tops\n")
        for name, path in sorted(tiles.items()):
            stem = os.path.splitext(name)[0]
            with rasterio.open(path) as src:
                raw = src.read(1).astype(np.float64)
                nod = src.nodata
            if nod is not None:
                raw[raw == nod] = 0.0
            raw[~np.isfinite(raw)] = 0.0
            sm = smooth_chm(raw, ws=SMOOTH_WS, method="median")
            sub = detect_tops(sm, hmin=HMIN, ws=WS)
            tops = as_pixels(sub)
            if not tops:
                continue
            g = HierarchicalRegionGrower(sm)
            out = g.run_all(tops, variance_thresh=VT, mask_thresh=MASK,
                            conflict_rule=RULE)
            pre = f"{OUT}/real_{stem}"
            np.savetxt(f"{pre}_chm.csv", raw, delimiter=",", fmt="%.17g")
            np.savetxt(f"{pre}_sm.csv", sm, delimiter=",", fmt="%.17g")
            np.savetxt(f"{pre}_sub.csv", np.asarray(sub, float) + 1.0,
                       delimiter=",", fmt="%.17g")
            np.savetxt(f"{pre}_mask.csv", g._prepare_mask(MASK).astype(int),
                       fmt="%d", delimiter=",")
            np.savetxt(f"{pre}_ws.csv", g.watershed_labels, fmt="%d", delimiter=",")
            np.savetxt(f"{pre}_out.csv", out, fmt="%d", delimiter=",")
            f.write(f"{stem},{len(tops)}\n")
    print(f"stage-by-stage reference for the real tiles written to {OUT}/real.csv")
