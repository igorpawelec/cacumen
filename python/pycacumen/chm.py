"""
pycacumen.chm — canopy height model smoothing.

A raw CHM carries pit noise and spurious local maxima: a single laser
return through a gap, or the ragged upper surface of a crown, both create
peaks that tree detection would read as separate trees. Smoothing trades
a little height accuracy for far fewer false tops.

Copyright (C) 2025 Igor Pawelec
Licence: GPLv3 — see LICENSE.
"""

import numpy as np
from numba import njit

__all__ = ["smooth_chm", "SMOOTHING_METHODS"]

SMOOTHING_METHODS = ("median", "mean", "gaussian", "maximum")


def _ensure_scipy():
    try:
        from scipy import ndimage
    except ImportError as e:  # pragma: no cover
        raise ImportError(
            "scipy is required for CHM smoothing.\n"
            "  conda install -c conda-forge scipy"
        ) from e
    return ndimage


def smooth_chm(chm, ws=3, method="median"):
    """
    Smooth a canopy height model.

    Parameters
    ----------
    chm : ndarray, shape (rows, cols)
        Canopy height model.
    ws : int
        Window size in pixels. Larger windows suppress more false tops but
        also merge genuinely adjacent crowns. Should be smaller than the
        narrowest crown you want to keep.
    method : {'median', 'mean', 'gaussian', 'maximum'}
        - ``'median'`` (default) — removes pits and spikes while keeping
          crown edges sharp. The usual choice.
        - ``'mean'`` — cheaper, but blurs crown boundaries.
        - ``'gaussian'`` — smooth falloff, sigma = ws / 3.
        - ``'maximum'`` — dilates crowns; flattens the apex, so tops are
          detected as plateaus rather than points. Rarely what you want
          before :func:`pycacumen.treetops.detect_tops`.

    Returns
    -------
    smoothed : ndarray, same shape as *chm*

    Notes
    -----
    Cells with no data (NaN) are skipped: a cell takes the statistic of the
    cells of its window that have data. The cells on the rim of a hole get a
    height that way, so a hole shrinks by ``ws // 2`` on each side, and a
    cell whose whole window has no data stays NaN.
    """
    ndimage = _ensure_scipy()
    chm = np.asarray(chm)
    if chm.ndim != 2:
        raise ValueError(f"chm must be 2-D, got shape {chm.shape}")
    if ws < 1:
        raise ValueError(f"ws must be >= 1, got {ws}")
    if method not in SMOOTHING_METHODS:
        raise ValueError(
            f"method must be one of {SMOOTHING_METHODS}, got {method!r}"
        )
    # An even window has no centre pixel, so scipy places it half a pixel
    # off, and the result then depends on which way the raster happens to
    # be oriented: smoothing a scene and smoothing its mirror image gave
    # values up to 8.8 m apart on a 40x55 test scene. Nobody picks ws=4
    # intending a shifted window, so this is refused rather than warned
    # about. 'gaussian' is exempt because ws only scales sigma there and
    # the kernel stays symmetric.
    if method != "gaussian" and ws % 2 == 0:
        raise ValueError(
            f"ws must be odd for method='{method}', got {ws}. An even window "
            f"has no centre pixel, so it sits half a pixel off and the "
            f"result depends on the raster's orientation. Use {ws - 1} or "
            f"{ws + 1}."
        )

    # NaN is skipped, not propagated and not fed to the filters raw.
    #
    # scipy's comparison-based filters are undefined on NaN rather than
    # NaN-aware: the result depends on where in the window the NaN sits.
    # With the nine values 1..9 and one NaN, median_filter returned 9, 5 or
    # 4 depending on the position, and maximum_filter returned nan, 9 or 8.
    # On chm_150_2023.tif that left 1208 pixels (3% of the raster) smoothed
    # to an arbitrary window element, 20 of 253 tree tops standing on them,
    # and errors against a NaN-skipping median of up to 26 m.
    #
    # The fast path is taken whenever there is no NaN, so a clean CHM costs
    # nothing and gives byte-identical results to before.
    if not np.isnan(chm).any():
        if method == "median":
            return ndimage.median_filter(chm, size=ws)
        if method == "mean":
            return ndimage.uniform_filter(chm, size=ws)
        if method == "gaussian":
            return ndimage.gaussian_filter(chm, sigma=ws / 3.0)
        return ndimage.maximum_filter(chm, size=ws)

    return _smooth_with_nan(ndimage, chm, ws, method)


# A window of up to this many values is sorted where it lies; of a larger one
# only the middle is found. Sorting by insertion costs the square of the count:
# nothing for a 5 x 5 window, and what made a 41 x 41 one slower in 0.7.1 than
# the code that 0.7.1 replaced.
_SORT_UP_TO = 25


@njit(cache=True, nogil=True)
def _kth_smallest(a, n, k):
    """Rearrange a[:n] so that a[k] is its k-th smallest value (from 0), with
    nothing larger before it and nothing smaller after it, and return a[k].

    Hoare's selection: the range that holds position k is split around one
    of its values until the split falls on k. The time grows with n, not
    with its square, also where all the values are equal, as on the ground.

    The value to split around is the one in the middle of the range for
    the first two rounds: the best guess where the values lie in order, as
    a slope read row by row does. Where they rise and then fall, as a
    ridge read row by row does, the middle holds the largest every time
    and each round takes off one row of the window, tens of times the work
    for 41 x 41 cells. From the third round the value comes from a place
    scattered over the range, which the orders a raster produces do not
    follow. The place is a function of the round alone, so the same window
    always takes the same course, and an order built against that very
    sequence can still be slow; it is not one a raster holds.
    """
    lo, hi = 0, n - 1
    rounds = 0
    step = 0
    while lo < hi:
        if rounds < 2:
            pivot = a[(lo + hi) // 2]
        else:
            step = (step * 1103515245 + 12345) & 0x7FFFFFFF
            pivot = a[lo + step % (hi - lo + 1)]
        rounds += 1
        i, j = lo, hi
        while i <= j:
            while a[i] < pivot:
                i += 1
            while a[j] > pivot:
                j -= 1
            if i <= j:
                a[i], a[j] = a[j], a[i]
                i += 1
                j -= 1
        if k <= j:
            hi = j
        elif k >= i:
            lo = i
        else:                                # between the two parts: a[k] is the pivot
            break
    return a[k]


@njit(cache=True, nogil=True)
def _nanmedian_filter(a, ws):
    """Median of the cells with data in a ws x ws window around every cell.

    The window is reflected at the edges as scipy.ndimage reflects it
    (d c b a | a b c d | d c b a), folding as often as it has to on a
    raster narrower than the window. An even count takes the mean of its
    two middle values, as numpy's median does. A window with no data gives
    NaN. Compiled, and without the interpreter lock.

    A window of a few values is sorted by insertion; of a larger one the
    middle value is selected and the rest left unsorted (_kth_smallest),
    so that the cost per cell grows with the window and not with its
    square.
    """
    rows, cols = a.shape
    half = ws // 2
    out = np.empty((rows, cols), dtype=np.float64)
    window = np.empty(ws * ws, dtype=np.float64)
    for r in range(rows):
        for c in range(cols):
            n = 0
            for dr in range(-half, half + 1):
                rr = (r + dr) % (2 * rows)
                if rr >= rows:
                    rr = 2 * rows - 1 - rr
                for dc in range(-half, half + 1):
                    cc = (c + dc) % (2 * cols)
                    if cc >= cols:
                        cc = 2 * cols - 1 - cc
                    v = a[rr, cc]
                    if not np.isnan(v):
                        window[n] = v
                        n += 1
            if n == 0:
                out[r, c] = np.nan
            elif n > _SORT_UP_TO:
                k = n // 2
                upper = _kth_smallest(window, n, k)
                if n % 2:
                    out[r, c] = upper
                else:                            # the other middle value: the largest before k
                    lower = window[0]
                    for i in range(1, k):
                        if window[i] > lower:
                            lower = window[i]
                    out[r, c] = (lower + upper) / 2.0
            else:
                for i in range(1, n):            # insertion sort: the window is small
                    v = window[i]
                    j = i - 1
                    while j >= 0 and window[j] > v:
                        window[j + 1] = window[j]
                        j -= 1
                    window[j + 1] = v
                if n % 2:
                    out[r, c] = window[n // 2]
                else:
                    out[r, c] = (window[n // 2 - 1] + window[n // 2]) / 2.0
    return out


def _smooth_with_nan(ndimage, chm, ws, method):
    """Window statistics over the non-NaN cells; all-NaN windows stay NaN."""
    valid = ~np.isnan(chm)
    filled = np.where(valid, chm, 0.0)
    w = valid.astype(np.float64)

    if method == "maximum":
        # -inf loses every comparison, so it drops out of the max without
        # poisoning it the way NaN does.
        out = ndimage.maximum_filter(np.where(valid, chm, -np.inf), size=ws)
    elif method == "mean":
        total = ndimage.uniform_filter(filled, size=ws)
        count = ndimage.uniform_filter(w, size=ws)
        with np.errstate(invalid="ignore", divide="ignore"):
            out = total / count
    elif method == "gaussian":
        # Normalised convolution: weight the kernel by which cells exist.
        sigma = ws / 3.0
        total = ndimage.gaussian_filter(filled, sigma=sigma)
        count = ndimage.gaussian_filter(w, sigma=sigma)
        with np.errstate(invalid="ignore", divide="ignore"):
            out = total / count
    else:                                    # median
        # A compiled pass: for every cell, the median of the cells of its
        # window that have data. Until 0.7.1 this was scipy's generic_filter
        # calling numpy's nanmedian once per cell, in Python: hundreds of
        # times slower than the median of a CHM without holes (a second or
        # more for a 200 x 201 tile), and holding the interpreter throughout,
        # so a program that smoothed in a worker thread stood still. The values
        # are the same, the cast back to the CHM's own precision included:
        # generic_filter stored its result in the type of its input. Only the
        # sign of a zero can differ, in a window that holds 0.0 and -0.0:
        # which of the two is the median was numpy's partition's to say.
        out = _nanmedian_filter(np.ascontiguousarray(chm, dtype=np.float64), int(ws))
        if chm.dtype.kind == "f" and chm.dtype.itemsize < 8:
            out = out.astype(chm.dtype).astype(np.float64)

    # A window with nothing in it has no statistic. For the gaussian that
    # means the *kernel* reached nothing, and the kernel is wider than ws --
    # 4 sigma either side -- so the test is the denominator rather than a
    # ws-sized window. Using ws there would have marked 7663 cells empty
    # against the 5889 the kernel actually fails to reach, and put rcacumen and
    # this package 1774 pixels apart.
    if method == "gaussian":
        empty = count == 0
    else:
        empty = ndimage.maximum_filter(valid.astype(np.uint8), size=ws) == 0
    out = np.asarray(out, dtype=float)
    out[empty] = np.nan
    out[~np.isfinite(out) & ~empty] = np.nan
    return out
