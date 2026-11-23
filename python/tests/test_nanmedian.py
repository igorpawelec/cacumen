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

from pycacumen.chm import smooth_chm, _nanmedian_filter

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
