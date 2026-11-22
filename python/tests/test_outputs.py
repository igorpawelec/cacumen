"""What the README's quick start writes, and what the command line prints.

Three things the quick start did not do as written until 0.6.2:
``cd.to_vector("out/")`` failed unless the folder already existed, the
``.vrt`` it wrote beside the polygons placed the raster above its own top
edge and read as zeros, and ``pycacumen -V`` asked for an input file.

Copyright (C) 2026 Igor Pawelec. Licence: GPLv3.
"""
import os
import subprocess
import sys
import warnings

import numpy as np
import pytest

rasterio = pytest.importorskip("rasterio")
fiona = pytest.importorskip("fiona")

from pycacumen import CrownDelineator, __version__, smooth_chm  # noqa: E402
from pycacumen.cli import build_parser  # noqa: E402

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, "test_data")
CHM = os.path.join(DATA, "chm_150_2014.tif")
needs_data = pytest.mark.skipif(not os.path.exists(CHM), reason="test_data/chm_150_2014.tif is not in this checkout")


@pytest.fixture(scope="module")
def delineated():
    cd = CrownDelineator.from_file(CHM, quiet=True)
    cd.smooth(ws=3).detect(hmin=7, ws=5).merge(5.0).screen(10.0)
    cd.delineate(variance_thresh=2.0)
    return cd


@needs_data
def test_to_vector_makes_the_folder_it_is_given(delineated, tmp_path):
    out = tmp_path / "out" / "nested"                       # neither exists
    delineated.to_vector(str(out) + "/", name="crowns")
    assert sorted(p.name for p in out.iterdir() if p.suffix in (".bin", ".vrt", ".shp")) == [
        "crowns.bin", "crowns.shp", "crowns.vrt"]
    with fiona.open(str(out / "crowns.shp")) as src:
        assert len(src) == len(set(np.unique(delineated.crowns)) - {0})


@needs_data
def test_tops_to_vector_makes_the_folder_too(delineated, tmp_path):
    out = tmp_path / "tops"
    delineated.tops_to_vector(str(out), name="tops")
    assert (out / "tops_treetops.shp").exists()


@needs_data
def test_the_vrt_is_the_crown_raster(delineated, tmp_path):
    """Same grid as the GeoTIFF, same values: it is the crowns as a raw file
    with a header, for software that reads no GeoTIFF."""
    delineated.to_raster(str(tmp_path / "crowns.tif"))
    delineated.to_vector(str(tmp_path / "v"), name="crowns")
    try:
        rasterio.open(str(tmp_path / "v" / "crowns.vrt")).close()
    except rasterio.errors.RasterioIOError as e:
        if "VRTRawRasterBand" in str(e) or "local file" in str(e):
            pytest.skip(f"this GDAL is configured not to read raw VRT sources: {e}")
        raise
    with rasterio.open(str(tmp_path / "crowns.tif")) as tif, rasterio.open(str(tmp_path / "v" / "crowns.vrt")) as vrt:
        assert vrt.bounds == tif.bounds
        assert vrt.transform == tif.transform
        assert vrt.bounds.bottom < vrt.bounds.top
        assert vrt.crs == tif.crs
        band = vrt.read(1)
        assert band.dtype == np.int32
        assert np.array_equal(band, tif.read(1))
        assert len(np.unique(band)) > 100                   # not the zeros it used to read as


def test_version_flag_needs_no_input(capsys):
    for flag in ("-V", "--version"):
        with pytest.raises(SystemExit) as stop:
            build_parser().parse_args([flag])
        assert stop.value.code == 0
        assert capsys.readouterr().out.strip() == f"pycacumen {__version__}"


def test_version_flag_from_the_command_line():
    r = subprocess.run([sys.executable, "-m", "pycacumen", "-V"], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert r.stdout.strip() == f"pycacumen {__version__}"


def test_smoothing_a_chm_with_holes_warns_about_nothing():
    """The all-NaN window is the expected case and NaN is its answer; numpy's
    "All-NaN slice encountered" for each of them was printed on every run.

    The warnings are recorded, not turned into errors: an exception raised
    inside the callback of scipy's generic_filter is swallowed there."""
    chm = np.full((30, 30), 20.0)
    chm[5:15, 5:15] = np.nan
    with warnings.catch_warnings(record=True) as seen:
        warnings.simplefilter("always")
        sm = smooth_chm(chm, ws=3, method="median")
    assert [str(w.message) for w in seen] == []
    assert np.isnan(sm[10, 10]) and sm[0, 0] == 20.0
    assert np.isnan(sm).sum() == 8 * 8                      # the hole shrinks by one cell at its rim


@needs_data
def test_max_height_is_the_highest_cell_under_the_crown(delineated, tmp_path):
    """The raw CHM has cells with no data inside the crowns, where smoothing
    filled the rim of a hole. np.maximum carries the NaN of one such cell to
    the whole crown: 124 of the 212 polygons of this tile had no max_height."""
    delineated.to_vector(str(tmp_path), name="crowns")
    with fiona.open(str(tmp_path / "crowns.shp")) as src:
        got = {f["properties"]["id"]: f["properties"]["max_height"] for f in src}
    crowns, chm = delineated.crowns, delineated.chm
    assert int((np.isnan(chm) & (crowns > 0)).sum()) > 0        # the tile does have such cells
    assert len(got) == 212
    wholly_on_holes = 0
    for crown_id, height in got.items():
        assert height is not None and np.isfinite(height), crown_id
        cells = chm[crowns == crown_id]
        if np.isnan(cells).all():
            # a crown of a few cells can stand on filled holes only: the smoothed surface is its height
            wholly_on_holes += 1
            cells = delineated.smoothed[crowns == crown_id]
        assert height == pytest.approx(float(np.nanmax(cells)), abs=0.006), crown_id
    assert wholly_on_holes == 2


@needs_data
def test_an_empty_folder_still_means_the_current_directory(delineated, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    delineated.to_vector("", name="here")
    assert (tmp_path / "here.shp").exists() and (tmp_path / "here.vrt").exists()


def test_a_crown_in_several_pieces_carries_the_area_of_each_piece(tmp_path):
    """One of the 145 crowns of this tile polygonises into three pieces. Each
    carried the area and the diameter of the whole crown, so the areas of the
    layer added up to more than the crowns cover."""
    tile = os.path.join(DATA, "chm_80_2021.tif")
    if not os.path.exists(tile):
        pytest.skip("test_data/chm_80_2021.tif is not in this checkout")
    cd = CrownDelineator.from_file(tile, quiet=True)
    cd.smooth(ws=3).detect(hmin=7, ws=5).merge(5.0).screen(10.0)
    crowns = cd.delineate(variance_thresh=2.0)
    cd.to_vector(str(tmp_path), name="crowns")
    with fiona.open(str(tmp_path / "crowns.shp")) as src:
        rows = [(f["properties"]["id"], f["properties"]["area_m2"], f["properties"]["crown_diam"]) for f in src]
    ids = [i for i, _, _ in rows]
    assert len(set(ids)) == len(set(np.unique(crowns)) - {0})
    assert len(ids) > len(set(ids))                          # the tile does have a crown in pieces
    pixel_area = abs(cd.transform.a * cd.transform.e)
    assert sum(a for _, a, _ in rows) == pytest.approx(int((crowns > 0).sum()) * pixel_area, abs=0.005 * len(rows))
    for _, area, diam in rows:
        assert diam == pytest.approx(2.0 * np.sqrt(area / np.pi), abs=0.011)
