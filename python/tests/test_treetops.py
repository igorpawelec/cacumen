"""Tests for pycacumen.chm and pycacumen.treetops."""
import numpy as np
import pytest

from pycacumen.chm import smooth_chm, SMOOTHING_METHODS
from pycacumen.treetops import detect_tops, merge_tops, screen_tops, as_pixels


class TestSmoothing:

    def test_all_methods_run(self, stand):
        chm, _ = stand
        for m in SMOOTHING_METHODS:
            out = smooth_chm(chm, ws=3, method=m)
            assert out.shape == chm.shape

    def test_median_reduces_spikes(self, stand):
        chm, _ = stand
        spiky = chm.copy()
        spiky[60, 60] = 200.0
        assert smooth_chm(spiky, ws=3, method="median")[60, 60] < 100.0

    def test_bad_method(self, stand):
        chm, _ = stand
        with pytest.raises(ValueError, match="method must be"):
            smooth_chm(chm, method="nope")

    def test_bad_ws(self, stand):
        chm, _ = stand
        with pytest.raises(ValueError, match="ws must be"):
            smooth_chm(chm, ws=0)

    def test_rejects_3d(self):
        with pytest.raises(ValueError, match="2-D"):
            smooth_chm(np.zeros((3, 10, 10)))


class TestDetectTops:

    def test_finds_the_trees(self, stand):
        chm, expected = stand
        tops = detect_tops(smooth_chm(chm, ws=3), hmin=5, ws=5)
        assert len(tops) == len(expected)

    def test_shape_and_dtype(self, stand):
        chm, _ = stand
        tops = detect_tops(chm, hmin=5, ws=5)
        assert isinstance(tops, np.ndarray)
        assert tops.ndim == 2 and tops.shape[1] == 2
        assert tops.dtype == np.float64

    def test_empty_is_still_2d(self, flat):
        """Regression: an empty (0,) array breaks any [:, 0] downstream."""
        tops = detect_tops(flat, hmin=20, ws=5)
        assert tops.shape == (0, 2)
        tops[:, 0]  # must not raise

    def test_hmin_filters(self, stand):
        chm, _ = stand
        assert len(detect_tops(chm, hmin=100, ws=5)) == 0

    def test_larger_window_detects_no_more(self, stand):
        chm, _ = stand
        sm = smooth_chm(chm, ws=3)
        assert len(detect_tops(sm, hmin=5, ws=9)) <= len(detect_tops(sm, hmin=5, ws=3))


class TestMergeTops:

    def test_merges_close_pair(self):
        tops = np.array([[10.0, 10.0], [10.0, 12.0], [50.0, 50.0]])
        merged = merge_tops(tops, distance=5.0)
        assert len(merged) == 2

    def test_keeps_distant(self):
        tops = np.array([[10.0, 10.0], [50.0, 50.0]])
        assert len(merge_tops(tops, distance=5.0)) == 2

    def test_merged_is_centroid(self):
        tops = np.array([[10.0, 10.0], [10.0, 20.0]])
        merged = merge_tops(tops, distance=20.0)
        assert len(merged) == 1
        np.testing.assert_allclose(merged[0], [10.0, 15.0])

    def test_empty_and_single(self):
        assert merge_tops(np.empty((0, 2))).shape == (0, 2)
        assert merge_tops(np.array([[1.0, 2.0]])).shape == (1, 2)

    def test_accepts_list(self):
        assert merge_tops([(1.0, 2.0), (1.0, 3.0)], distance=5.0).shape == (1, 2)


class TestScreenTops:

    def test_drops_short(self, stand):
        chm, _ = stand
        sm = smooth_chm(chm, ws=3)
        tops = detect_tops(sm, hmin=5, ws=5)
        kept = screen_tops(sm, tops, hmin=100.0)
        assert len(kept) == 0
        assert kept.shape == (0, 2)

    def test_keeps_tall(self, stand):
        chm, _ = stand
        sm = smooth_chm(chm, ws=3)
        tops = detect_tops(sm, hmin=5, ws=5)
        assert len(screen_tops(sm, tops, hmin=0.0)) == len(tops)

    def test_out_of_bounds_top_is_dropped_not_wrapped(self, stand):
        """Negative indices would silently wrap around the raster."""
        chm, _ = stand
        kept = screen_tops(chm, np.array([[-5.0, -5.0], [500.0, 500.0]]), hmin=0.0)
        assert len(kept) == 0


class TestAsPixels:

    def test_a_top_goes_to_the_nearest_cell(self):
        """Whole numbers are the centres of the cells. Until 0.7.0 the
        coordinate was cut down, and (3.7, 4.2) went to cell (3, 4)."""
        px = as_pixels(np.array([[3.7, 4.2]]))
        assert px == [(4, 4)]
        assert all(isinstance(v, int) for v in px[0])
        assert as_pixels([[186.67, 48.67], [0.49, 0.0], [2.0, 7.0]]) == [(187, 49), (0, 0), (2, 7)]
        assert as_pixels([[0.5, 1.5], [2.5, 3.5]]) == [(1, 2), (3, 4)]     # a half goes up

    def test_the_last_bit_of_a_plateau_centre_does_not_choose_the_cell(self):
        """A plateau two cells long has its centre on the boundary between them. The
        weighted mean that finds it gives 40.5 or 40.49999999999999 according to the
        order of the sums, and without a margin rcacumen and pycacumen seeded
        different cells on a double-precision CHM."""
        for height, first in ((12.3, 40), (7.7, 6), (0.1 + 0.2, 10)):
            chm = np.zeros((60, 20))
            chm[first:first + 2, 4] = height
            tops = detect_tops(chm, hmin=0.05, ws=3)
            assert len(tops) == 1 and abs(tops[0, 0] - (first + 0.5)) < 1e-9
            assert as_pixels(tops) == [(first + 1, 4)]
        assert as_pixels([[40.5 - 1e-9, 2.0], [40.5 - 1e-3, 2.0]]) == [(41, 2), (40, 2)]

    def test_empty(self):
        assert as_pixels(np.empty((0, 2))) == []


class TestEvenWindowRejected:
    """An even window has no centre pixel and sits half a pixel off.

    Measured before the guard went in: smoothing a 40x55 scene and its
    mirror image differed by up to 8.8 m at ws=4, and detect_tops found
    397 tops against 400 on the mirror of chm_150_2023.tif at ws=4, and
    188 against 206 at ws=6. Nobody chooses ws=4 intending a shifted
    window, and nothing in the output shows that it happened.
    """

    @staticmethod
    def _scene():
        rng = np.random.default_rng(4)
        return rng.random((40, 55)) * 20 + 5

    @pytest.mark.parametrize("method", ["median", "mean", "maximum"])
    @pytest.mark.parametrize("ws", [2, 4, 6])
    def test_smooth_rejects_even_ws(self, method, ws):
        with pytest.raises(ValueError, match="odd"):
            smooth_chm(self._scene(), ws=ws, method=method)

    def test_smooth_allows_even_ws_for_gaussian(self):
        """ws only scales sigma there, so the kernel stays symmetric."""
        sc = self._scene()
        a = smooth_chm(sc, ws=4, method="gaussian")
        b = smooth_chm(sc[:, ::-1], ws=4, method="gaussian")[:, ::-1]
        np.testing.assert_allclose(a, b)

    @pytest.mark.parametrize("ws", [2, 4, 6])
    def test_detect_rejects_even_ws(self, ws):
        with pytest.raises(ValueError, match="odd"):
            detect_tops(self._scene(), hmin=5, ws=ws)

    @pytest.mark.parametrize("method", ["median", "mean", "maximum", "gaussian"])
    def test_odd_ws_is_orientation_independent(self, method):
        sc = self._scene()
        a = smooth_chm(sc, ws=5, method=method)
        b = smooth_chm(sc[:, ::-1], ws=5, method=method)[:, ::-1]
        np.testing.assert_allclose(a, b)

    def test_detection_with_odd_ws_is_orientation_independent(self):
        sc = self._scene()
        assert len(detect_tops(sc, hmin=5, ws=5)) == \
            len(detect_tops(sc[:, ::-1], hmin=5, ws=5))


class TestOrder:
    """rcacumen holds itself to the same two orders (r/tests/test-chm-treetops.R),
    which is what makes a crown id mean the same tree in both packages."""

    def test_tops_come_in_the_order_of_a_row_by_row_reading(self):
        chm = np.zeros((30, 30))
        chm[4, 19] = 20.0
        chm[9, 4] = 25.0
        assert detect_tops(chm, hmin=5, ws=3).tolist() == [[4.0, 19.0], [9.0, 4.0]]

    def test_merged_groups_come_in_the_order_of_their_first_member(self):
        tops = np.array([[10, 10], [30, 30], [11, 10], [50, 50], [31, 30]], float)
        assert merge_tops(tops, distance=2.0).tolist() == [[10.5, 10.0], [30.5, 30.0], [50.0, 50.0]]
