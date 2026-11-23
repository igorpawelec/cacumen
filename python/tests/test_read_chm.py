"""A raster's nodata value is read as NaN, as rcacumen reads it as NA."""
import numpy as np
import pytest

from pycacumen import CrownDelineator


class TestReadChm:

    rasterio = pytest.importorskip("rasterio")

    def _write(self, path, data, nodata, dtype):
        from rasterio.transform import from_origin
        with self.rasterio.open(str(path), "w", driver="GTiff", height=data.shape[0],
                                width=data.shape[1], count=1, dtype=dtype, nodata=nodata,
                                crs="EPSG:2180", transform=from_origin(500000, 600000, 0.5, 0.5)) as dst:
            dst.write(data.astype(dtype), 1)

    @pytest.mark.parametrize("dtype,nodata", [("float32", -9999.0), ("int16", -32768)])
    def test_a_declared_nodata_value_is_read_as_nan(self, tmp_path, dtype, nodata):
        from pycacumen.io import read_chm
        data = np.full((12, 15), 20.0)
        data[3:6, 4:9] = nodata
        self._write(tmp_path / "chm.tif", data, nodata, dtype)
        chm, transform, crs = read_chm(tmp_path / "chm.tif")
        assert chm.dtype == np.float32 and chm.shape == (12, 15)
        assert np.isnan(chm[3:6, 4:9]).all() and np.isnan(chm).sum() == 15
        assert np.nanmin(chm) == 20.0
        win, _, _ = read_chm(tmp_path / "chm.tif", window=(4, 3, 6, 4))
        assert win.shape == (4, 6) and np.isnan(win[0:3, 0:5]).all() and win[3, 0] == 20.0

    def test_nan_as_nodata_and_no_nodata_at_all_read_as_before(self, tmp_path):
        from pycacumen.io import read_chm
        data = np.full((8, 8), 12.5)
        data[2, 2] = np.nan
        self._write(tmp_path / "nan.tif", data, float("nan"), "float32")
        chm, _, _ = read_chm(tmp_path / "nan.tif")
        assert np.isnan(chm[2, 2]) and np.isnan(chm).sum() == 1
        data[2, 2] = -9999.0
        self._write(tmp_path / "none.tif", data, None, "float32")
        chm, _, _ = read_chm(tmp_path / "none.tif")
        assert chm[2, 2] == -9999.0 and not np.isnan(chm).any()   # nothing declares it

    def test_from_file_smooths_past_the_hole_as_rcacumen_does(self, tmp_path):
        """Read as a height, -9999 entered the median of the cells beside the hole."""
        data = np.full((9, 9), 20.0)
        data[:, 4:] = 24.0
        data[3:6, 3:6] = -9999.0
        self._write(tmp_path / "hole.tif", data, -9999.0, "float32")
        cd = CrownDelineator.from_file(tmp_path / "hole.tif", quiet=True).smooth(ws=3)
        assert np.nanmin(cd.smoothed) >= 20.0
        assert np.isnan(cd.smoothed[4, 4])

    def test_a_mask_band_is_not_consulted(self, tmp_path):
        """rcacumen reads through terra, which goes by the nodata value alone;
        the two packages have to see the same holes."""
        from pycacumen.io import read_chm
        from rasterio.transform import from_origin
        data = np.full((6, 6), 15.0, dtype="float32")
        mask = np.full((6, 6), 255, dtype="uint8")
        mask[1:3, 1:3] = 0
        with self.rasterio.Env(GDAL_TIFF_INTERNAL_MASK=True):
            with self.rasterio.open(str(tmp_path / "masked.tif"), "w", driver="GTiff", height=6,
                                    width=6, count=1, dtype="float32", crs="EPSG:2180",
                                    transform=from_origin(0, 6, 1, 1)) as dst:
                dst.write(data, 1)
                dst.write_mask(mask)
        chm, _, _ = read_chm(tmp_path / "masked.tif")
        assert not np.isnan(chm).any() and (chm == 15.0).all()
