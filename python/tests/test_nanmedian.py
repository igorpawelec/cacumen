"""The median of a CHM with holes, compiled.

Until 0.7.1 a CHM with no data in it was smoothed by scipy's
``generic_filter`` calling numpy's ``nanmedian`` once per cell: correct,
hundreds of times slower than the median of a CHM without holes, and
holding the interpreter for the whole of it. The compiled pass that took
its place has to return the same numbers, on every shape and in every
precision, because the tops and the crowns are found on them. Equal
numbers are equal to the last bit, with one exception: 0.0 and -0.0.

Copyright (C) 2026 Igor Pawelec. Licence: GPLv3.
"""
import warnings

import numpy as np
import pytest

import pycacumen.chm as chm_module
from pycacumen.chm import smooth_chm, _kth_smallest, _nanmedian_filter, _SORT_UP_TO

ndimage = pytest.importorskip("scipy.ndimage")


def _reference(chm, ws):
    """What 0.7.0 computed: nanmedian per window, stored in the type of the input."""
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="All-NaN slice", category=RuntimeWarning)
        out = ndimage.generic_filter(chm, np.nanmedian, size=ws, mode="reflect")
    return out.astype(np.float64)


def _with_holes(rng, shape, dtype, share):
    chm = rng.uniform(0.0, 40.0, shape).astype(dtype)
    chm[rng.random(shape) < share] = np.nan
    if not np.isnan(chm).any():                  # the compiled pass runs only on a CHM with holes
        chm.flat[0] = np.nan
    return chm


class TestSameAsNanmedianPerWindow:

    @pytest.mark.parametrize("dtype", [np.float32, np.float64])
    @pytest.mark.parametrize("ws", [1, 3, 5, 7])
    def test_on_rasters_with_holes(self, dtype, ws):
        rng = np.random.default_rng(ws)
        for shape in ((23, 31), (40, 17), (8, 8)):
            for share in (0.05, 0.4, 0.9):
                chm = _with_holes(rng, shape, dtype, share)
                got = smooth_chm(chm, ws=ws, method="median")
                assert got.dtype == np.float64
                assert np.array_equal(got, _reference(chm, ws), equal_nan=True), (shape, share)

    @pytest.mark.parametrize("dtype", [np.float32, np.float64])
    @pytest.mark.parametrize("ws", [7, 11, 17, 21])
    def test_on_windows_too_large_to_sort(self, dtype, ws):
        """Above 25 values the middle of the window is selected, not sorted for.
        Whole metres make ties, of which a real CHM has many: the ground."""
        assert ws * ws > _SORT_UP_TO
        rng = np.random.default_rng(100 + ws)
        for shape in ((40, 37), (25, 60)):
            for share in (0.05, 0.5):
                for whole in (False, True):
                    chm = _with_holes(rng, shape, dtype, share)
                    if whole:
                        chm = np.round(chm / 8.0)
                    got = smooth_chm(chm, ws=ws, method="median")
                    assert np.array_equal(got, _reference(chm, ws), equal_nan=True), \
                        (shape, share, whole)

    def test_a_window_many_times_wider_than_the_raster(self):
        """The reflected edge folds four times over where the window is eight times
        the raster. scipy reads outside the array there, which is what 0.7.0
        returned on a CHM with holes, so the reference is written out here."""
        def folded(chm, ws):
            rows, cols = chm.shape
            half = ws // 2
            rr = np.arange(-half, rows + half) % (2 * rows)
            cc = np.arange(-half, cols + half) % (2 * cols)
            padded = chm[np.where(rr < rows, rr, 2 * rows - 1 - rr)][
                :, np.where(cc < cols, cc, 2 * cols - 1 - cc)]
            out = np.full(chm.shape, np.nan)
            for r in range(rows):
                for c in range(cols):
                    cells = padded[r:r + ws, c:c + ws]
                    if not np.isnan(cells).all():
                        out[r, c] = np.nanmedian(cells)
            return out

        rng = np.random.default_rng(23)
        for shape, ws in (((2, 21), 17), ((28, 2), 19), ((4, 9), 33), ((3, 3), 25)):
            chm = _with_holes(rng, shape, np.float64, 0.3)
            got = smooth_chm(chm, ws=ws, method="median")
            assert np.array_equal(got, folded(chm, ws), equal_nan=True), (shape, ws)

    @pytest.mark.parametrize("shape", [(1, 13), (13, 1), (2, 3), (1, 1), (3, 2)])
    def test_on_rasters_narrower_than_the_window(self, shape):
        """The reflected edge folds more than once there."""
        rng = np.random.default_rng(7)
        for ws in (3, 5, 9):
            chm = _with_holes(rng, shape, np.float64, 0.3)
            got = smooth_chm(chm, ws=ws, method="median")
            assert np.array_equal(got, _reference(chm, ws), equal_nan=True), (shape, ws)

    def test_single_precision_keeps_its_own_values(self):
        """An even count is the mean of two values, which a float32 cannot always
        hold; generic_filter stored it in the type of the CHM, and so do we."""
        chm = np.array([[0.1, np.nan, 0.2], [np.nan, np.nan, np.nan], [0.4, np.nan, 0.7]],
                       dtype=np.float32)
        got = smooth_chm(chm, ws=3, method="median")
        assert np.array_equal(got, got.astype(np.float32).astype(np.float64))
        assert np.array_equal(got, _reference(chm, 3), equal_nan=True)


    def test_zeros_of_both_signs_give_a_zero(self):
        """Which of 0.0 and -0.0 is the median of a window that holds both was up
        to numpy's partition, and is up to the order of the cells now. The two
        compare equal and nothing downstream tells them apart; the sign is the
        one thing that is not held to 0.7.0."""
        chm = np.array([[-0.0, 0.0, 5.0], [np.nan, 1.0, -0.0], [0.0, -0.0, 2.0]])
        got = smooth_chm(chm, ws=3, method="median")
        assert np.array_equal(got, _reference(chm, 3), equal_nan=True)
        assert got[0, 0] == 0.0


class TestTheMedianItself:

    def test_an_even_count_takes_the_mean_of_the_middle_two(self):
        chm = np.array([[1.0, np.nan, 5.0], [np.nan, 2.0, np.nan], [8.0, np.nan, np.nan]])
        assert smooth_chm(chm, ws=3, method="median")[1, 1] == 3.5      # of 1, 2, 5, 8

    def test_a_window_with_no_data_stays_empty_and_says_nothing(self):
        chm = np.full((9, 9), 20.0)
        chm[2:9, 2:9] = np.nan
        with warnings.catch_warnings():
            warnings.simplefilter("error")                              # numpy warned per cell
            got = smooth_chm(chm, ws=3, method="median")
            nothing = smooth_chm(np.full((4, 5), np.nan), ws=3, method="median")
        assert np.isnan(got[5, 5]) and got[0, 0] == 20.0
        assert np.isnan(got[4:8, 4:8]).all() and np.isfinite(got[2, 2])  # the rim is filled
        assert nothing.shape == (4, 5) and np.isnan(nothing).all()

    def test_the_input_is_left_as_it_was(self):
        chm = _with_holes(np.random.default_rng(1), (12, 12), np.float32, 0.3)
        before = chm.copy()
        smooth_chm(chm, ws=5, method="median")
        assert np.array_equal(chm, before, equal_nan=True)

    def test_it_runs_without_the_interpreter_lock(self):
        """A program that smooths in a worker thread keeps its window alive only
        if the pass lets the interpreter go."""
        assert _nanmedian_filter.targetoptions.get("nogil") is True
        assert _kth_smallest.targetoptions.get("nogil") is True


class TestSelection:
    """The middle of a large window is found without sorting the window."""

    def test_the_kth_value_stands_at_k_with_the_smaller_before_it(self):
        rng = np.random.default_rng(11)
        for n in list(range(1, 40)) + [81, 225, 441, 1681]:
            for kind in ("spread", "ties", "equal", "rising", "falling"):
                values = {"spread": rng.uniform(0, 40, n),
                          "ties": rng.integers(0, 4, n).astype(float),
                          "equal": np.full(n, 7.5),
                          "rising": np.arange(n, dtype=float),
                          "falling": np.arange(n, dtype=float)[::-1].copy()}[kind]
                ordered = np.sort(values)
                for k in {0, n // 2, max(0, n // 2 - 1), n - 1, int(rng.integers(0, n))}:
                    a = np.concatenate([values, [-1.0, 99.0]])       # beyond n: not to be touched
                    assert _kth_smallest(a, n, k) == ordered[k], (n, kind, k)
                    assert a[k] == ordered[k]
                    assert (a[:k] <= a[k]).all() and (a[k + 1:n] >= a[k]).all(), (n, kind, k)
                    assert np.array_equal(np.sort(a[:n]), ordered)   # the same values, rearranged
                    assert a[n] == -1.0 and a[n + 1] == 99.0

    def test_windows_of_25_and_of_26_values(self):
        """25 values are sorted and 26 selected from: the two ways meet here. The
        window of the middle cell of a 7 x 7 raster is the raster, so the number
        of values in it is the number of cells left with data."""
        assert _SORT_UP_TO == 25
        rng = np.random.default_rng(3)
        for count in (23, 24, 25, 26, 27, 28, 49):
            for whole in (False, True):
                for _ in range(20):
                    chm = rng.uniform(0.0, 40.0, (7, 7))
                    if whole:
                        chm = np.round(chm / 8.0)
                    empty = rng.permutation(49)[:49 - count]
                    chm.flat[empty] = np.nan
                    values = chm[~np.isnan(chm)]
                    assert len(values) == count
                    got = _nanmedian_filter(chm, 7)[3, 3]
                    assert got == np.median(values), (count, whole)

    def test_selection_is_used_above_25_values_and_only_there(self, monkeypatch):
        """The uncompiled source of the pass, run with a selection that counts its
        calls: one call for every window of more than 25 values, none for the
        others, and the values of the compiled pass."""
        rng = np.random.default_rng(8)
        chm = rng.uniform(0.0, 40.0, (9, 8))
        chm[rng.random(chm.shape) < 0.45] = np.nan

        ws, half = 7, 3                          # the number of values in every window
        rr = np.arange(-half, 9 + half) % 18
        cc = np.arange(-half, 8 + half) % 16
        data = ~np.isnan(chm)[np.where(rr < 9, rr, 17 - rr)][:, np.where(cc < 8, cc, 15 - cc)]
        counts = np.array([[data[r:r + ws, c:c + ws].sum() for c in range(8)] for r in range(9)])
        assert (counts > _SORT_UP_TO).any() and (counts <= _SORT_UP_TO).any()

        seen = []

        def counting(a, n, k):
            seen.append(n)
            return _kth_smallest.py_func(a, n, k)

        monkeypatch.setattr(chm_module, "_kth_smallest", counting)
        got = _nanmedian_filter.py_func(chm, ws)
        assert sorted(seen) == sorted(counts[counts > _SORT_UP_TO].tolist())
        monkeypatch.undo()
        assert np.array_equal(got, _nanmedian_filter(chm, ws), equal_nan=True)

    def test_the_orders_a_raster_produces_do_not_make_the_selection_slow(self):
        """Split around the middle value every time, a window whose values rise and
        then fall, a ridge read row by row, loses one row a round: half a million
        reads for 41 x 41 values, where twelve thousand do. Counted on the
        uncompiled source, so that no clock is involved."""
        class Counted:
            def __init__(self, values):
                self.values = list(values)
                self.reads = 0

            def __getitem__(self, i):
                self.reads += 1
                return self.values[i]

            def __setitem__(self, i, value):
                self.values[i] = value

        n = 41 * 41
        half = n // 2
        yy, xx = np.mgrid[0:41, 0:41]
        against_the_middle = np.empty(n)         # the smallest left always in the middle
        places = list(range(n))
        for value in range(n):
            against_the_middle[places.pop((len(places) - 1) // 2)] = value
        orders = {
            "random": np.random.default_rng(1).uniform(0, 40, n),
            "rising": np.arange(n, dtype=float),
            "falling": np.arange(n, dtype=float)[::-1],
            "rising then falling": np.concatenate([np.arange(half, dtype=float),
                                                  np.arange(n - half, dtype=float)[::-1]]),
            "all equal": np.zeros(n),
            "a ridge, row by row": (-np.abs(yy - 20) + 0.001 * xx).ravel(),
            "a dome, row by row": (30.0 - 0.05 * ((yy - 18.6) ** 2 + (xx - 24.1) ** 2)).ravel(),
            "built against the middle": against_the_middle,
        }
        for name, values in orders.items():
            a = Counted(values)
            assert _kth_smallest.py_func(a, n, half) == np.sort(values)[half], name
            assert a.reads < 12 * n, (name, a.reads)
