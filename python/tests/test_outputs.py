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
from pycacumen.treetops import _cells  # noqa: E402

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
@pytest.mark.parametrize("driver, ext", [("ESRI Shapefile", ".shp"), ("GPKG", ".gpkg")])
def test_a_top_carries_the_height_it_was_detected_at(delineated, tmp_path, driver, ext):
    """Until 0.7.1 ``height`` was read from the CHM as given. A top is found on
    the smoothed one, and the cell under it is often a hole or a pit that the
    smoothing filled: 22 of the 215 tops of the quick start had no height or a
    height of 0.00 m, one of them 12.98 m high where it was kept \"from 10 m\".
    ``height`` is the smoothed value now and ``raw_height`` the one that was
    there before."""
    delineated.tops_to_vector(str(tmp_path), name="tops", driver=driver)
    with fiona.open(str(tmp_path / ("tops_treetops" + ext))) as src:
        assert list(src.schema["properties"]) == ["id", "height", "raw_height"]
        rows = [dict(f["properties"]) for f in src]
    assert [r["id"] for r in rows] == list(range(215))

    cells = _cells(delineated.tops)
    smoothed = delineated.smoothed[cells[:, 0], cells[:, 1]]
    raw = delineated.chm[cells[:, 0], cells[:, 1]]
    for row, s, r in zip(rows, smoothed, raw):
        assert row["height"] == pytest.approx(float(s), abs=0.006)
        if np.isnan(r):
            assert row["raw_height"] is None
        else:
            assert row["raw_height"] == pytest.approx(float(r), abs=0.006)

    assert min(r["height"] for r in rows) > 10.0             # what screen(10.0) kept
    assert sum(r["raw_height"] is None for r in rows) == 21
    pits = [r for r in rows if r["raw_height"] is not None and r["raw_height"] <= 10.0]
    assert [(r["raw_height"], r["height"]) for r in pits] == [(0.0, 12.98)]


def test_tops_written_by_hand_carry_what_they_are_given(tmp_path):
    """``save_tree_tops`` on its own writes ``height`` from the CHM it is handed
    and no second field, as it always did; a cell with no data is null, where
    it used to be written as not-a-number."""
    from affine import Affine
    from pycacumen.io import save_tree_tops

    chm = np.arange(20.0).reshape(4, 5)
    chm[2, 3] = np.nan
    tops = np.array([[1.0, 1.0], [2.0, 3.0], [0.6, 3.4]])
    wkt = rasterio.crs.CRS.from_epsg(2180).to_wkt()
    transform = Affine(0.5, 0.0, 500000.0, 0.0, -0.5, 300000.0)

    save_tree_tops(tops, str(tmp_path), "plain", transform, wkt, chm, driver="GPKG")
    with fiona.open(str(tmp_path / "plain_treetops.gpkg")) as src:
        assert list(src.schema["properties"]) == ["id", "height"]
        assert [f["properties"]["height"] for f in src] == [6.0, None, 8.0]

    filled = np.where(np.isnan(chm), 12.5, chm)
    save_tree_tops(tops, str(tmp_path), "both", transform, wkt, filled, driver="GPKG",
                   raw_chm=chm)
    with fiona.open(str(tmp_path / "both_treetops.gpkg")) as src:
        got = [(f["properties"]["height"], f["properties"]["raw_height"]) for f in src]
    assert got == [(6.0, 6.0), (12.5, None), (8.0, 8.0)]


def test_tops_put_in_by_hand_on_a_chm_that_was_not_smoothed_carry_its_own_height(tmp_path):
    """detect() smooths a CHM that was not smoothed, so there is a smoothed CHM
    whenever the tops came from the pipeline. Without one the CHM as given is
    all there is."""
    from affine import Affine

    yy, xx = np.mgrid[0:30, 0:30]
    chm = (20.0 * np.exp(-((yy - 15) ** 2 + (xx - 15) ** 2) / 40.0)).astype(np.float32)
    cd = CrownDelineator(chm, transform=Affine(1.0, 0.0, 0.0, 0.0, -1.0, 30.0),
                         crs=rasterio.crs.CRS.from_epsg(2180), quiet=True)
    cd.tops = np.array([[15.0, 15.0]])
    assert cd.smoothed is None
    cd.tops_to_vector(str(tmp_path), name="t", driver="GPKG")
    with fiona.open(str(tmp_path / "t_treetops.gpkg")) as src:
        got = [dict(f["properties"]) for f in src]
    assert got == [{"id": 0, "height": 20.0, "raw_height": 20.0}]


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

    The warnings are recorded, not turned into errors: while the median was
    scipy's generic_filter, an exception raised inside its callback was
    swallowed there."""
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
    assert len(got) == len(np.unique(crowns[crowns > 0]))
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
