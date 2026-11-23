"""
pycacumen.treetops — tree top detection, merging and screening.

Tree tops are found as local maxima of a smoothed CHM, the standard
raster approach to individual tree detection (Popescu & Wynne 2004). The
window size sets the implicit minimum spacing between trees, so it is the
main lever on over- versus under-detection.

Every function here returns tree tops as an ``(n, 2)`` float64 array of
(row, column) coordinates — including the empty case, which is ``(0, 2)``
rather than ``(0,)`` so that callers can index ``[:, 0]`` unconditionally.

References
----------
Popescu, S.C., Wynne, R.H. (2004). Seeing the trees in the forest: using
    lidar and multispectral data fusion with local filtering and variable
    window size for estimating tree height. Photogrammetric Engineering &
    Remote Sensing 70(5), 589-604.

Copyright (C) 2025 Igor Pawelec
Licence: GPLv3 — see LICENSE.
"""

import warnings

import numpy as np

__all__ = ["detect_tops", "merge_tops", "screen_tops", "as_pixels"]


def _ensure_scipy():
    try:
        from scipy import ndimage
    except ImportError as e:  # pragma: no cover
        raise ImportError(
            "scipy is required for tree top detection.\n"
            "  conda install -c conda-forge scipy"
        ) from e
    return ndimage


def _as_tops_array(x):
    """Normalise anything top-like to float64 (n, 2), empty -> (0, 2)."""
    return np.asarray(x, dtype=np.float64).reshape(-1, 2)


def detect_tops(chm, hmin=2.0, ws=3):
    """
    Detect tree tops as local maxima of a (smoothed) CHM.

    A pixel is a candidate when it equals the maximum of its ``ws x ws``
    neighbourhood and exceeds *hmin*. Adjacent candidates — a flat apex
    yields several — are labelled into connected components and each
    component is reduced to its height-weighted centre of mass, giving
    subpixel coordinates.

    Parameters
    ----------
    chm : ndarray, shape (rows, cols)
        Canopy height model. Smooth it first (:func:`pycacumen.chm.smooth_chm`);
        on a raw CHM this will over-detect badly.
    hmin : float
        Minimum height for a pixel to be considered. Sets the floor
        between canopy and understorey/ground.
    ws : int
        Neighbourhood size in pixels. Acts as the minimum spacing between
        detected tops: too small over-detects one tree as many, too large
        merges neighbouring trees. Over-detection is the safer error here,
        since :func:`pycacumen.hrg.HierarchicalRegionGrower.run_all` can merge
        surplus tops back together.

    Returns
    -------
    tops : ndarray, shape (n, 2), dtype float64
        Subpixel (row, column) coordinates. Empty gives shape (0, 2). They
        come in the order of the first cell of each maximum when the raster
        is read row by row; rcacumen lists them in the same order.
    """
    ndimage = _ensure_scipy()
    chm = np.asarray(chm)
    if chm.ndim != 2:
        raise ValueError(f"chm must be 2-D, got shape {chm.shape}")
    if ws < 1:
        raise ValueError(f"ws must be >= 1, got {ws}")
    # An even window has no centre pixel, so the "local maximum" test is
    # taken over a window sitting half a pixel off, and the tree count then
    # depends on the raster's orientation: on chm_150_2023.tif, ws=4 found
    # 397 tops and 400 on the mirror image, ws=6 found 188 and 206. Refused
    # rather than warned about, because the failure is invisible in the
    # output.
    if ws % 2 == 0:
        raise ValueError(
            f"ws must be odd, got {ws}. An even window has no centre pixel, "
            f"so the local-maximum test sits half a pixel off and the tree "
            f"count depends on the raster's orientation. Use {ws - 1} or "
            f"{ws + 1}."
        )

    # -inf, not NaN, for the cells with no data. scipy's maximum_filter is
    # undefined on NaN rather than NaN-aware: with the values 1..9 and one
    # NaN it returned nan, 9 or 8 depending on where the NaN sat. -inf loses
    # every comparison and drops out cleanly.
    #
    # A NaN pixel is never itself detected: NaN == anything is False, and so
    # is NaN > hmin.
    local_max = ndimage.maximum_filter(
        np.where(np.isnan(chm), -np.inf, chm), size=ws)
    detected = (chm == local_max) & (chm > hmin)
    labels, num = ndimage.label(detected)
    centres = ndimage.center_of_mass(chm, labels, range(1, num + 1))
    return _as_tops_array(centres)


# A merged group that stretches over more than this many merge distances was
# joined as a chain, tree to tree, and merge_tops says so.
CHAIN_FACTOR = 5.0

# How close to the boundary between two cells a top counts as on it (see _cells).
_HALF_MARGIN = 1e-6


def _cells(tops):
    """The cell each top stands on, as an (n, 2) integer array.

    A top is a position in array coordinates, where whole numbers are the
    centres of the cells, so the cell under it is the nearest one: 186.67
    stands in cell 187. Until 0.7.0 the coordinate was cut down, to 186,
    which put the seed of a plateau one cell up or left of its centre for
    nearly half of the tops, and now and then on the ground beside the
    crown.

    A half goes up, the same way in rcacumen, and so does anything within
    a millionth of a pixel of a half. The centre of a plateau two cells
    long is a half exactly, but it is computed as a weighted mean, and
    comes out as 40.5 or as 40.49999999999999 according to the order of
    the sums: without the margin the last bit chose the cell, and chose
    differently in the two packages.
    """
    return np.floor(np.asarray(tops, dtype=float) + (0.5 + _HALF_MARGIN)).astype(np.intp)


def _heights_at(chm, tops):
    """CHM value under each top; NaN for a top outside the raster or on no data."""
    rows, cols = chm.shape
    rc = _cells(tops)
    inside = ((rc[:, 0] >= 0) & (rc[:, 0] < rows) &
              (rc[:, 1] >= 0) & (rc[:, 1] < cols))
    h = np.full(len(tops), np.nan)
    h[inside] = chm[rc[inside, 0], rc[inside, 1]]
    return h


def merge_tops(tops, distance=5.0, chm=None, hmin=None):
    """
    Merge tree tops that lie within *distance* of each other.

    Tops closer than the threshold are almost always the same tree seen
    through several local maxima. Neighbours within the radius are found
    with a k-d tree, grouped transitively with union-find, and each group
    is replaced by its centroid.

    A centroid is a point between the tops. On a CHM with holes it can fall
    on a cell with no data, and between two crowns it can fall in the gap,
    on the ground. Nothing grows from there, or a crown grows on the
    ground: the trees the group stood for are lost. Given the CHM, such a
    group is represented by the highest of its own tops instead.

    A warning is given when groups stretch over more than five times the
    merge distance: those were joined as chains, tree to tree, and the
    distance is too large for the stand.

    Note the grouping is transitive: a chain of tops each within
    *distance* of the next collapses into one, even if the endpoints are
    far apart. Keep the threshold well below the typical crown diameter.

    Parameters
    ----------
    tops : array-like, shape (n, 2)
        Tree top coordinates.
    distance : float
        Merge radius in pixels.
    chm : ndarray, shape (rows, cols), optional
        The CHM the tops were detected on. With it, a group whose centroid
        falls on a cell with no data (NaN) or outside the raster is replaced
        by the highest of its tops (the first of them when several are
        equally high). Without it every group is replaced by its centroid,
        wherever that falls.
    hmin : float, optional
        The height the tops were detected above. With it and the CHM, a
        group whose centroid stands on a cell that is not above *hmin*,
        where no top could have been detected, is replaced by its highest
        top as well.

    Returns
    -------
    merged : ndarray, shape (m, 2), m <= n
        The groups in the order of their first member in *tops*; rcacumen
        returns them in the same order.
    """
    try:
        from scipy.spatial import cKDTree
    except ImportError as e:  # pragma: no cover
        raise ImportError(
            "scipy is required for merging tree tops.\n"
            "  conda install -c conda-forge scipy"
        ) from e

    if not (np.isfinite(distance) and distance >= 0):
        raise ValueError(f"distance must be a number, zero or more, got {distance!r}")
    if hmin is not None and not np.isfinite(hmin):
        raise ValueError(f"hmin must be a number, got {hmin!r}")
    tops = _as_tops_array(tops)
    if tops.shape[0] < 2:
        return tops

    kd = cKDTree(tops)
    neighbour_lists = kd.query_ball_point(tops, r=distance)

    n = tops.shape[0]
    parent = np.arange(n)

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]  # path compression
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    for i, neighbours in enumerate(neighbour_lists):
        for j in neighbours:
            if j > i:
                union(i, j)

    from collections import defaultdict
    groups = defaultdict(list)
    for i in range(n):
        groups[find(i)].append(i)

    members = list(groups.values())
    merged = _as_tops_array([tops[idx].mean(axis=0) for idx in members])

    # Grouping is transitive, so tops a merge distance apart, one after
    # another, become one group however long the row. With a distance near
    # the spacing of the trees a whole stand folds into a few tops, and
    # nothing in the result says that it happened.
    spans = np.array([np.hypot(*np.ptp(tops[idx], axis=0)) for idx in members])
    far = spans > CHAIN_FACTOR * distance
    if far.any():
        worst = int(np.argmax(spans))
        warnings.warn(
            f"merge_tops: {int(far.sum())} of {len(members)} groups stretch over more than "
            f"{CHAIN_FACTOR:g} times the merge distance; the longest joins "
            f"{len(members[worst])} tops over {spans[worst]:.1f} px. Tops merge in chains, "
            f"each within {distance:g} px of the next, so a distance near the spacing of "
            "the trees runs neighbouring trees together. Lower the distance.",
            stacklevel=2)

    if chm is not None:
        chm = np.asarray(chm)
        h_tops = _heights_at(chm, tops)
        h_merged = _heights_at(chm, merged)
        lost = np.isnan(h_merged)
        if hmin is not None:
            # not above hmin: where detect_tops could not have put a top
            with np.errstate(invalid="ignore"):
                lost |= ~(h_merged > hmin)
        for g in np.flatnonzero(lost):
            h = h_tops[members[g]]
            if len(h) > 1 and not np.isnan(h).all():
                merged[g] = tops[members[g][int(np.nanargmax(h))]]
    return merged


def screen_tops(chm, tops, hmin):
    """
    Drop tree tops standing on CHM values below *hmin*.

    Useful as a second pass with a stricter threshold than
    :func:`detect_tops` used — e.g. detect at 2 m to let the growing see
    the understorey, then keep only tops above 10 m.

    Parameters
    ----------
    chm : ndarray, shape (rows, cols)
        The same CHM the tops were detected on.
    tops : array-like, shape (n, 2)
    hmin : float
        Minimum height to keep.

    Returns
    -------
    kept : ndarray, shape (m, 2), m <= n
    """
    chm = np.asarray(chm)
    tops = _as_tops_array(tops)
    if tops.shape[0] == 0:
        return tops

    rows, cols = chm.shape
    rc = _cells(tops)
    inside = ((rc[:, 0] >= 0) & (rc[:, 0] < rows) &
              (rc[:, 1] >= 0) & (rc[:, 1] < cols))
    keep = np.zeros(len(tops), dtype=bool)
    keep[inside] = chm[rc[inside, 0], rc[inside, 1]] >= hmin
    return tops[keep]


def as_pixels(tops):
    """
    Convert subpixel tops to integer (row, col) pixel indices.

    Each top goes to the cell it stands on, the nearest one: whole numbers
    are the centres of the cells, and a half goes up.

    Returns
    -------
    pixels : list of (int, int)
        Suitable for :meth:`pycacumen.hrg.HierarchicalRegionGrower.run_all`.
    """
    tops = _as_tops_array(tops)
    return [(int(r), int(c)) for r, c in _cells(tops)]
