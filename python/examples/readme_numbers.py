"""
readme_numbers.py
~~~~~~~~~~~~~~~~~
The numbers the README quotes in "Notes on behaviour" and "Performance",
measured on synthetic stands that this script builds, so that they can be
measured again.

Usage, from anywhere::

    python examples/readme_numbers.py

A stand is a square CHM with one Gaussian crown per cell of a jittered grid,
random in height and width, smoothed and searched for tops with the package's
own functions. The package is read from the checkout this script sits in.
Timings are the best of a few repeats, after one call that compiles the
kernels; they are this machine's and are quoted in the README as round
figures. Nothing is written.
"""
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

from pycacumen import (smooth_chm, detect_tops, as_pixels,          # noqa: E402
                       HierarchicalRegionGrower)
from pycacumen.hrg import _build_adjacency_and_stats               # noqa: E402


def stand(size, spacing, seed):
    """A size x size CHM with a crown in every cell of a jittered grid."""
    rng = np.random.default_rng(seed)
    yy, xx = np.mgrid[0:size, 0:size]
    chm = np.zeros((size, size))
    for r0 in range(spacing // 2, size, spacing):
        for c0 in range(spacing // 2, size, spacing):
            r = r0 + rng.uniform(-spacing / 4, spacing / 4)
            c = c0 + rng.uniform(-spacing / 4, spacing / 4)
            h = rng.uniform(12, 30)
            s = rng.uniform(spacing / 5, spacing / 3)
            win = (slice(max(0, int(r - 4 * s)), int(r + 4 * s) + 1),
                   slice(max(0, int(c - 4 * s)), int(c + 4 * s) + 1))
            chm[win] = np.maximum(chm[win], h * np.exp(
                -((yy[win] - r) ** 2 + (xx[win] - c) ** 2) / (2 * s ** 2)))
    sm = smooth_chm(chm.astype(np.float32), ws=3, method="median")
    return sm, as_pixels(detect_tops(sm, hmin=5.0, ws=5))


def best(f, repeats=3):
    times = []
    for _ in range(repeats):
        t0 = time.perf_counter()
        out = f()
        times.append(time.perf_counter() - t0)
    return min(times), out


def grow(sm, tops, **kw):
    return HierarchicalRegionGrower(sm).run_all(tops, mask_thresh=1.0, **kw)


def main():
    # compile the kernels on a small stand, outside every timing
    sm0, tops0 = stand(60, 12, 0)
    grow(sm0, tops0, variance_thresh=8.0)

    print("retry_rejected (README, 'only bites in the middle'):")
    for seed in (1, 2, 3):
        sm, tops = stand(500, 15, seed)
        line = []
        for vt in (2.0, 8.0, 20.0, 120.0):
            a = grow(sm, tops, variance_thresh=vt, retry_rejected=False)
            b = grow(sm, tops, variance_thresh=vt, retry_rejected=True)
            na, nb = len(np.unique(a)) - 1, len(np.unique(b)) - 1
            line.append("vt %g: %.1f%% of the raster, crowns %d -> %d"
                        % (vt, 100.0 * (a != b).mean(), na, nb))
        print("  stand %d, %d x %d px, %d tops: %s"
              % (seed, sm.shape[0], sm.shape[1], len(tops), "; ".join(line)))

    print("growing (README, 'n_jobs' and 'Performance'):")
    sm, tops = stand(500, 15, 1)
    t1, seq = best(lambda: grow(sm, tops, variance_thresh=8.0))
    t2, par = best(lambda: grow(sm, tops, variance_thresh=8.0, n_jobs=2))
    print("  %d tops on %d x %d px: run_all %.2f s sequential, %.2f s with n_jobs=2; "
          "results identical: %s" % (len(tops), sm.shape[0], sm.shape[1], t1, t2,
                                     bool(np.array_equal(seq, par))))

    print("the pixel pass (README, 'Performance'):")
    sm, tops = stand(1000, 14, 1)
    g = HierarchicalRegionGrower(sm)
    markers = np.zeros(sm.shape, dtype=np.int32)
    for i, (r, c) in enumerate(tops, start=1):
        markers[r, c] = i
    g.initial_watershed(markers, mask=g._prepare_mask(1.0))
    tk, _ = best(lambda: _build_adjacency_and_stats(g.watershed_labels, g.chm,
                                                    g.num_regions), repeats=5)
    tg, _ = best(lambda: g.build_adjacency_graph(), repeats=3)
    tw, _ = best(lambda: grow(sm, tops, variance_thresh=8.0), repeats=2)
    print("  %d tops on %d x %d px: the compiled pass over the pixels %.1f ms, "
          "the whole graph %.2f s, the whole run_all %.1f s"
          % (len(tops), sm.shape[0], sm.shape[1], 1000.0 * tk, tg, tw))

if __name__ == "__main__":
    main()
