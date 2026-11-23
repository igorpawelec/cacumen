"""A merged top that lands where no tree stands, on a hole or in the gap between two
crowns, does not take its trees with it; and a merge that runs through the stand in
chains says so."""
import warnings

import numpy as np
import pytest

from pycacumen import CrownDelineator
from pycacumen.treetops import merge_tops, screen_tops, _cells


def _two_crowns(n=41, left=17, right=23, hole=True):
    """Two crowns on row 20, the right-hand one the taller. With ``hole``, a block
    of no data between them, wide enough to outlast the smoothing, which fills
    the cells on a hole's rim."""
    yy, xx = np.mgrid[0:n, 0:n]
    chm = np.zeros((n, n))
    for r, c, h in ((20, left, 22.0), (20, right, 26.0)):
        chm = np.maximum(chm, h * np.exp(-((yy - r) ** 2 + (xx - c) ** 2) / (2 * 2.5 ** 2)))
    if hole:
        chm[18:23, 19:22] = np.nan
    return chm.astype(np.float32)


class TestMergedTopOnAHole:

    TOPS = np.array([[20.0, 17.0], [20.0, 23.0]])          # whole numbers: the centres of two cells

    def test_without_the_chm_the_centroid_stands_wherever_it_falls(self):
        merged = merge_tops(self.TOPS, distance=7.0)
        assert merged.shape == (1, 2) and np.allclose(merged[0], [20.0, 20.0])

    def test_with_the_chm_the_highest_top_stands_in(self):
        chm = _two_crowns()
        assert np.isnan(chm[20, 20])                        # where the centroid falls
        merged = merge_tops(self.TOPS, distance=7.0, chm=chm)
        assert merged.shape == (1, 2) and np.allclose(merged[0], [20.0, 23.0])
        # and it survives the screening that used to drop the centroid
        assert len(screen_tops(chm, merged, hmin=10.0)) == 1
        assert len(screen_tops(chm, merge_tops(self.TOPS, distance=7.0), hmin=10.0)) == 0

    def test_equally_high_tops_give_the_first(self):
        chm = _two_crowns()
        chm[20, 17] = chm[20, 23] = 30.0
        merged = merge_tops(self.TOPS, distance=7.0, chm=chm)
        assert np.allclose(merged[0], [20.0, 17.0])
        merged = merge_tops(self.TOPS[::-1], distance=7.0, chm=chm)
        assert np.allclose(merged[0], [20.0, 23.0])

    def test_a_centroid_on_data_is_left_alone(self):
        chm = _two_crowns()
        chm[18:23, 20] = 5.0                                # the gap has data
        merged = merge_tops(self.TOPS, distance=7.0, chm=chm)
        assert np.allclose(merged[0], [20.0, 20.0])

    def test_groups_that_are_not_affected_keep_their_place_and_order(self):
        chm = _two_crowns()
        tops = np.array([[5.0, 5.0], [20.0, 17.0], [6.0, 5.0], [20.0, 23.0], [35.0, 35.0]])
        plain = merge_tops(tops, distance=7.0)
        aware = merge_tops(tops, distance=7.0, chm=chm)
        assert plain.shape == aware.shape == (3, 2)
        moved = np.any(plain != aware, axis=1)
        assert moved.sum() == 1 and np.allclose(aware[moved][0], [20.0, 23.0])

    def test_a_centroid_outside_the_raster_counts_as_no_data(self):
        chm = np.full((10, 10), 20.0, dtype=np.float32)
        tops = np.array([[9.0, 2.0], [12.0, 2.0]])          # the second is off the raster
        merged = merge_tops(tops, distance=5.0, chm=chm)
        assert np.allclose(merged[0], [9.0, 2.0])

    def test_the_pipeline_keeps_a_crown_for_the_group(self):
        chm = _two_crowns()
        cd = CrownDelineator(chm, quiet=True)
        cd.smooth(ws=3).detect(hmin=5, ws=3)
        assert len(cd.tops) == 2
        cd.merge(7.0).screen(10.0)
        assert len(cd.tops) == 1
        r, c = _cells(cd.tops)[0]
        assert np.isfinite(cd.smoothed[r, c]) and c > 20     # the taller, right-hand tree
        crowns = cd.delineate(variance_thresh=50.0, mask_thresh=1.0)
        assert (crowns > 0).sum() > 20


class TestMergedTopInAGap:
    """Two crowns far enough apart to leave ground between them: the centroid of their
    tops stands on that ground, where no top could have been detected."""

    TOPS = np.array([[20.0, 12.0], [20.0, 28.0]])

    def _chm(self):
        return _two_crowns(left=12, right=28, hole=False)

    def test_the_gap_is_ground(self):
        assert self._chm()[20, 20] < 0.5

    def test_without_hmin_the_centroid_stays_on_the_ground(self):
        merged = merge_tops(self.TOPS, distance=17.0, chm=self._chm())
        assert np.allclose(merged[0], [20.0, 20.0])

    def test_with_hmin_the_highest_top_stands_in(self):
        merged = merge_tops(self.TOPS, distance=17.0, chm=self._chm(), hmin=5.0)
        assert np.allclose(merged[0], [20.0, 28.0])

    def test_a_centroid_above_hmin_is_left_alone(self):
        chm = _two_crowns(left=18, right=22, hole=False)     # the crowns overlap: a saddle, not a gap
        tops = np.array([[20.0, 18.0], [20.0, 22.0]])
        assert chm[20, 20] > 5.0
        merged = merge_tops(tops, distance=5.0, chm=chm, hmin=5.0)
        assert np.allclose(merged[0], [20.0, 20.0])

    def test_hmin_is_the_threshold_of_detection_itself(self):
        """detect_tops takes cells above hmin; a centroid exactly at hmin is not above it."""
        chm = self._chm()
        chm[20, 20] = 5.0
        assert np.allclose(merge_tops(self.TOPS, distance=17.0, chm=chm, hmin=5.0)[0], [20.0, 28.0])
        chm[20, 20] = 5.5
        assert np.allclose(merge_tops(self.TOPS, distance=17.0, chm=chm, hmin=5.0)[0], [20.0, 20.0])

    def test_the_pipeline_passes_the_threshold_it_detected_with(self):
        cd = CrownDelineator(self._chm(), quiet=True)
        cd.smooth(ws=3).detect(hmin=5, ws=3)
        assert len(cd.tops) == 2
        cd.merge(17.0)
        assert len(cd.tops) == 1
        r, c = _cells(cd.tops)[0]
        assert cd.smoothed[r, c] > 5.0 and c > 20            # on the taller tree, not on the ground


class TestArguments:

    TOPS = np.array([[1.0, 1.0], [1.0, 3.0]])

    def test_a_distance_that_means_nothing_is_refused(self):
        for bad in (-2.0, float("nan"), float("inf")):
            with pytest.raises(ValueError, match="distance must be a number, zero or more"):
                merge_tops(self.TOPS, distance=bad)
        assert len(merge_tops(self.TOPS, distance=0.0)) == 2          # zero merges nothing

    def test_an_hmin_that_is_not_a_number_is_refused(self):
        with pytest.raises(ValueError, match="hmin must be a number"):
            merge_tops(self.TOPS, distance=3.0, chm=np.ones((5, 5)), hmin=float("nan"))


class TestChains:

    def test_a_merge_that_runs_in_a_chain_says_so(self):
        row = np.array([[10.0, 4.0 * i] for i in range(12)])           # each 4 px from the next
        with pytest.warns(UserWarning, match="stretch over more than 5 times the merge distance"):
            merged = merge_tops(row, distance=5.0)
        assert len(merged) == 1                                         # the result is what it was

    def test_the_warning_counts_the_groups_and_names_the_longest(self):
        long_row = [[10.0, 4.0 * i] for i in range(12)]                 # 44 px
        short_row = [[60.0, 4.0 * i] for i in range(8)]                 # 28 px
        pair = [[90.0, 10.0], [90.0, 13.0]]
        with pytest.warns(UserWarning) as record:
            merge_tops(np.array(long_row + short_row + pair), distance=5.0)
        message = str(record[0].message)
        assert "2 of 3 groups" in message and "joins 12 tops over 44.0 px" in message

    def test_compact_groups_are_merged_in_silence(self):
        tops = np.array([[10.0, 10.0], [12.0, 11.0], [11.0, 13.0], [40.0, 40.0], [42.0, 41.0]])
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            merged = merge_tops(tops, distance=5.0)
        assert len(merged) == 2
