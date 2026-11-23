"""A merged top that lands on a hole does not take its trees with it."""
import numpy as np

from pycacumen import CrownDelineator
from pycacumen.treetops import merge_tops, screen_tops


def _two_crowns_and_a_hole(n=41):
    """Two crowns six cells apart with a hole between them, wide enough to
    outlast the smoothing, which fills the cells on a hole's rim."""
    yy, xx = np.mgrid[0:n, 0:n]
    chm = np.zeros((n, n))
    for r, c, h in ((20, 17, 22.0), (20, 23, 26.0)):
        chm = np.maximum(chm, h * np.exp(-((yy - r) ** 2 + (xx - c) ** 2) / (2 * 2.5 ** 2)))
    chm[18:23, 19:22] = np.nan
    return chm.astype(np.float32)


class TestMergedTopOnAHole:

    TOPS = np.array([[20.5, 17.5], [20.5, 23.5]])

    def test_without_the_chm_the_centroid_stands_wherever_it_falls(self):
        merged = merge_tops(self.TOPS, distance=7.0)
        assert merged.shape == (1, 2) and np.allclose(merged[0], [20.5, 20.5])

    def test_with_the_chm_the_highest_top_stands_in(self):
        chm = _two_crowns_and_a_hole()
        assert np.isnan(chm[20, 20])                        # where the centroid falls
        merged = merge_tops(self.TOPS, distance=7.0, chm=chm)
        assert merged.shape == (1, 2) and np.allclose(merged[0], [20.5, 23.5])
        # and it survives the screening that used to drop the centroid
        assert len(screen_tops(chm, merged, hmin=10.0)) == 1
        assert len(screen_tops(chm, merge_tops(self.TOPS, distance=7.0), hmin=10.0)) == 0

    def test_equally_high_tops_give_the_first(self):
        chm = _two_crowns_and_a_hole()
        chm[20, 17] = chm[20, 23] = 30.0
        merged = merge_tops(self.TOPS, distance=7.0, chm=chm)
        assert np.allclose(merged[0], [20.5, 17.5])
        merged = merge_tops(self.TOPS[::-1], distance=7.0, chm=chm)
        assert np.allclose(merged[0], [20.5, 23.5])

    def test_a_centroid_on_data_is_left_alone(self):
        chm = _two_crowns_and_a_hole()
        chm[18:23, 20] = 5.0                                # the gap has data
        merged = merge_tops(self.TOPS, distance=7.0, chm=chm)
        assert np.allclose(merged[0], [20.5, 20.5])

    def test_groups_that_are_not_affected_keep_their_place_and_order(self):
        chm = _two_crowns_and_a_hole()
        tops = np.array([[5.5, 5.5], [20.5, 17.5], [6.5, 5.5], [20.5, 23.5], [35.5, 35.5]])
        plain = merge_tops(tops, distance=7.0)
        aware = merge_tops(tops, distance=7.0, chm=chm)
        assert plain.shape == aware.shape == (3, 2)
        moved = np.any(plain != aware, axis=1)
        assert moved.sum() == 1 and np.allclose(aware[moved][0], [20.5, 23.5])

    def test_a_centroid_outside_the_raster_counts_as_no_data(self):
        chm = np.full((10, 10), 20.0, dtype=np.float32)
        tops = np.array([[9.5, 2.5], [12.5, 2.5]])          # the second is off the raster
        merged = merge_tops(tops, distance=5.0, chm=chm)
        assert np.allclose(merged[0], [9.5, 2.5])

    def test_the_pipeline_keeps_a_crown_for_the_group(self):
        chm = _two_crowns_and_a_hole()
        cd = CrownDelineator(chm, quiet=True)
        cd.smooth(ws=3).detect(hmin=5, ws=3)
        assert len(cd.tops) == 2
        cd.merge(7.0).screen(10.0)
        assert len(cd.tops) == 1
        r, c = np.floor(cd.tops[0]).astype(int)
        assert np.isfinite(cd.smoothed[r, c]) and c > 20     # the taller, right-hand tree
        crowns = cd.delineate(variance_thresh=50.0, mask_thresh=1.0)
        assert (crowns > 0).sum() > 20
