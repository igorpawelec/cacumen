# Changelog


## [0.7.2] — 2026-11-23

### Fixed
- **The median of a CHM with holes no longer falls behind with a large window.** The
  compiled pass of 0.7.1 sorted every window by insertion, at a cost that grows with
  the square of the number of cells. With the windows of a few pixels that the method
  is used with, that is nothing. From a window of about 25 px the pass was slower than
  the code it had replaced, several times slower at 41 px (measured against the tag
  `python-v0.7.1`). A window of more than 25 values now has its middle value selected,
  and the rest is left unsorted. At 21 px the pass is several times ahead of numpy's
  `nanmedian` cell by cell, at 41 px less than twice (`examples/readme_numbers.py`).
  The values are those of 0.7.1 and of 0.7.0; as before, only the sign of a zero can
  differ.
- **A window many times wider than the raster.** Where the window is at least eight
  times as wide as the raster, a raster 2 cells wide under a window of 17 px for one,
  scipy reads outside the array when it reflects the edge, and on a CHM with holes
  0.7.0 returned what it found there. Since 0.7.1 the edge folds as often as it has
  to, which a test now holds (`tests/test_nanmedian.py`); the entry of 0.7.1 did not
  say so. A CHM without holes still goes through scipy.


## [0.7.1] — 2026-11-23

### Changed
- **A CHM with holes is smoothed by a compiled median.** Where the CHM had cells with
  no data, the median was scipy's `generic_filter` calling numpy's `nanmedian` once per
  cell: a second or two for 200 × 200 px with a quarter of the cells empty and tens of
  seconds for 1,000 × 1,000 px. It held the interpreter throughout, so a program that
  smoothed in a worker thread stood still for that long. The pass is compiled now and
  lets the interpreter go: milliseconds and a fraction of a second on the same rasters,
  over a hundred times less (`examples/readme_numbers.py`). The values are the same,
  in single and in double precision and on rasters narrower than the window
  (`tests/test_nanmedian.py`), so tops and crowns are what they were. Only the sign of
  a zero can differ, in a window that holds both 0.0 and -0.0.
- **The height beside an exported top is the height it was detected at.**
  `tops_to_vector` wrote `height` from the CHM as given, while a top is found, merged
  and screened on the smoothed one, and the cell under it is often a hole or a pit that
  the smoothing filled. With the quick start on `chm_150_2014.tif`, 21 of the 215 tops
  had no height and one had 0.00 m where the smoothed CHM has 12.98 m, though all of
  them were kept "from 10 m". `height` is read from the smoothed CHM now, and so meets
  the thresholds the top passed. The value of the CHM as given is written beside it as
  `raw_height`, null where the cell has no data. `save_tree_tops` takes `raw_chm` for
  that field, and writes a height that is not a number as null. Heights are rounded to
  two decimals in double precision now: 12.35 is written where rounding in the single
  precision of the CHM left 12.340000152587890. The positions of the tops are
  unchanged.
- **The chain warning of `merge_tops` is written for someone who has not read the
  code.** It speaks of chains of tops and no longer of groups. Every number comes with
  its unit: chains "longer than 25 px", where "more than 5 times the merge distance"
  stood beside a distance of "5 px". It points to the spacing of the detected tops and
  to the detection window, where it named the spacing of the trees: tops are detected
  closer together than the trees stand, the more so the smaller the window. And it has
  a sentence of its own for the case when all the tops were merged into one. The text
  is the same in rcacumen 0.6.1.
  When the warning is given has not changed.


## [0.7.0] — 2026-11-23

### Changed
- **A top goes to the nearest cell.** A top is a position in array coordinates, where
  whole numbers are the centres of the cells, and the cell under it was found by
  cutting the coordinate down: 186.67 went to cell 186, though it stands in 187. For a
  top that is the centre of a plateau, nearly half of the tops on the tiles of
  `test_data`, that put the seed one cell up or left of the centre, and now and then on
  the ground beside the crown, where a crown then grew: on `chm_11_2013.tif` the top at
  (186.67, 48.67) was seeded on a cell 0.05 m high next to one of 18.9 m. `as_pixels`,
  `screen_tops`, the height written beside an exported top and the fall-back of
  `merge_tops` use the nearest cell now. A half goes up, and so does anything within
  a millionth of a pixel of one: the centre of a plateau two cells long is a half that
  comes out with a last bit either way, and that bit must not choose the cell. On a
  plateau that bends, the nearest cell can still lie beside it. Tops and crowns change
  wherever a seed moves: the quick start on `chm_150_2014.tif` gives 215 tops and 214
  crowns where 0.6.3 gave 213 and 212.
  With the seeds on their plateaus the two packages also agree more closely. All
  seventy synthetic scenes are identical with rcacumen 0.6.0, where the ten large ones
  differed on 0.032 % of pixels, and on the seven real tiles the finished crowns differ
  on 0.10 % of pixels, where they differed on 0.64 %. The crowns are steadier as well:
  a nanometre of noise on a tile moves up to 24 % of its canopy pixels to another
  crown, where it moved up to 85 % (`r/tools/measure_tie_sensitivity.py`).
- **A merged group whose centroid stands where no top could be is represented by its
  highest top**, and not only a group whose centroid has no data. `merge_tops` takes
  `hmin`, the height the tops were detected above. With it and the CHM, a centroid on a
  cell that is not above `hmin`, the ground in the gap between two crowns for one, gives
  way to the highest top of the group. `CrownDelineator.merge` passes the threshold that
  `detect` used. On `chm_150_2014.tif` with `hmin` 2, a 3 px window and a merge distance
  of 5 px, three merged tops stood on cells 0.0 m high. `merge_tops` refuses a
  `distance` below zero and an `hmin` that is not a number.

### Added
- **`merge_tops` warns when it merges in chains.** Grouping is transitive, and with a
  distance near the spacing of the trees a whole stand folds into a few tops: 1,264
  tops into one at 8 px on `chm_33_2012.tif`. A `UserWarning` says how many groups
  stretch over more than five times the merge distance and how long the longest is. The
  result is unchanged.


## [0.6.3] — 2026-11-23

### Changed
- **`read_chm` and `CrownDelineator.from_file` return cells with no data as NaN.**
  The raster's nodata value used to come back as it was stored, so a declared
  -9999 was read as a height. It stayed out of every crown, being below any
  mask threshold, but it took part in the smoothing of the cells beside a hole,
  which rcacumen, reading `NA`, never let it do. Cells that hold the nodata
  value are NaN now, and smoothing and detection skip them. A raster whose
  nodata is NaN already, as the tiles in `test_data` are, reads as before, and
  so does one that declares none. A mask band or an alpha band is not
  consulted, as it is not in rcacumen.

### Fixed
- **A group of merged tops whose centroid fell on a hole left no crown.**
  `merge_tops` replaces tops closer than the radius by their centroid, and on a
  CHM with holes the centroid can fall on a cell with no data. The merged top
  was then screened out, or stood outside the mask, and the trees it stood for
  got no crown at all. `merge_tops` takes the CHM as an optional `chm` argument
  now and represents such a group by the highest of its own tops;
  `CrownDelineator.merge` and `delineate_crowns` pass it. On the tiles of
  `test_data` with the quick-start parameters one tree comes back on each of
  `chm_150_2023.tif` (214 to 215 tops) and `chm_33_2012.tif` (391 to 392); the
  other five are unchanged. `merge_tops(tops, distance)` without the CHM does
  what it did.


## [0.6.2] — 2026-11-22

### Fixed
- **`to_vector` and `--vector` failed unless the folder already existed.**
  The README's quick start writes to `out/`, and that raised
  `FileNotFoundError` on `out/crowns.bin`; from the command line the crown
  raster had been written by then and the run ended with exit code 1. The
  folder is created, for the tree-top points too.
- **The `.vrt` written beside the polygons did not describe the raster.** The
  pixel height was written with its sign reversed, which placed the raster
  above its own top edge (bottom 188,718.0 over top 188,617.5 on the test
  tile), and the band was declared as a plain VRT band, which has no source
  and reads as zeros. It is a `VRTRawRasterBand` on the `.bin` now, with the
  geotransform as it is, and opens as the same grid and values as the
  GeoTIFF.
- **`max_height` of the crown polygons was missing for most crowns of a CHM
  with holes.** A crown reaches over the rim of a hole, because smoothing
  gives those cells a height, and the maximum was taken with a function that
  propagates the NaN the raw CHM has there: with the quick-start parameters
  124 of the 212 polygons of
  `chm_150_2014.tif` had no height, 323 of 386 on `chm_33_2012.tif`. Cells
  with no data are left out of the maximum. A crown that stands on such
  cells only takes the highest cell of the smoothed CHM, the surface it was
  grown on. **This changes the `max_height` attribute**, not the crowns.
- **A crown that polygonises into several pieces carried the whole crown's
  area and diameter on each piece.** One of the 145 crowns of
  `chm_80_2021.tif` comes out in three, and the areas of the layer added up
  to more than the crowns cover. `area_m2` and `crown_diameter` are those of
  each polygon now, measured on its geometry.
- **`pycacumen -V` asked for `-i` and `-o`.** It prints the version and exits.
- Smoothing a CHM with holes no longer prints numpy's "All-NaN slice
  encountered" for every empty window; the result is unchanged.
- The tree-top points no longer use `Affine * (x, y)`, which affine 3
  announces the end of. Same coordinates.


## [0.6.1] — 2026-11-06

### Fixed
- **`run_all` raised `IndexError` when the last tree top had no watershed
  region.** A seed on a cell the mask removes (no data, below `mask_thresh`,
  cleaned away by `morpho_radius`) gets no region, so it has nothing to grow
  from. As the last seed its id lay past the region arrays and the run
  crashed; elsewhere in the list it grew an empty region. Such a seed now
  yields no crown, like an absorbed tree, its crown id stays unused, and
  `n_without_region` says how many there were. Found through OBIAStudio,
  where tree tops come from a CHM with no-data patches.


## [0.6.0] — 2026-08-24

### Changed
- **Renamed from `pyHRG`.** The import name is now `pycacumen`, so this is a
  breaking change for anything that imported the old one. Nothing else moved:
  the algorithms, the outputs and the cross-language agreement with the R twin
  are unchanged.


## [0.5.2] — 2026-07-22

### Fixed
- **A tree-top-out-of-bounds error crashed instead of explaining itself on a
  cp1250 console.** The message read `({rows}×{cols})` with `×` as U+00D7;
  on a Polish Windows console that raised `UnicodeEncodeError` while trying
  to raise the `ValueError`, replacing a clear diagnostic with a `charmap`
  traceback. The `×` is now `x`. `test_console_ascii.py` statically rejects
  non-ASCII in any `print`/`raise`/`warn` string, so the next one fails in CI
  rather than on a user's terminal.


## [0.5.1] — 2026-07-22

### Fixed
- **0.5.0 did not import at all on Python 3.9.** `max_iters: int | None` is
  PEP 604, which before 3.10 is a runtime expression: `TypeError: unsupported
  operand type(s) for |: 'type' and 'NoneType'`, raised while the class body
  executes, so `import pycacumen.hrg` failed outright. `pyproject.toml` claims
  `requires-python = ">=3.9"`. Now `Optional[int]`.

  It reached a tag because development runs on 3.12 and every local test
  passed. `tests/test_python39_compat.py` now checks the package the way 3.9
  would: no PEP 604 in annotations outside a module with
  `from __future__ import annotations`, and every file `compile()`s — the
  first attempt at this fix put the future import after three `__author__`
  assignments, where `ast.parse` accepts it and Python does not.


## [0.5.0] — 2026-07-22

### Fixed
- **NaN was handed straight to `scipy.ndimage`, which is undefined on it
  rather than NaN-aware.** With the nine values 1..9 and a single NaN,
  `median_filter` returned 9, 5 or 4 depending on *where* in the window the
  NaN sat, and `maximum_filter` returned nan, 9 or 8. A NaN-skipping median
  would have been a constant 5.5. The comparison-based selection is poisoned
  by NaN and picks an arbitrary element.

  On `chm_150_2023.tif` that left **1208 pixels — 3 % of the raster —
  smoothed to an arbitrary window value**, with **20 of 253 tree tops
  standing on them** and errors against a NaN-skipping median of up to
  **26.29 m**. Nothing raised, nothing warned; the output looked like a
  canopy.

  `smooth_chm()` and `detect_tops()` now skip NaN: window statistics are
  taken over the cells that exist, and a window holding nothing but NaN
  stays NaN. The fast path still runs whenever the CHM has no NaN, so a
  clean raster costs nothing and gives byte-identical results to before.

  Found while fixing rcacumen, which crashed on the same input. Neither was
  right: R failed loudly with an unhelpful message, Python produced
  plausible garbage. Both now agree, and the shared cross-check gained ten
  scenes with nodata holes so the path is actually covered.


## [0.4.0] — 2026-07-22

### Fixed
- **Exported tree tops sat half a pixel up and to the left.**
  `save_tree_tops()` mapped a top through `transform * (col, row)`, but an
  affine transform treats whole numbers as pixel *corners* while a top from
  `center_of_mass` is an array index, and array indices refer to pixel
  centres. Every exported point was therefore offset by half a pixel — 0.25 m
  on the 0.5 m test rasters, systematic, and enough to matter when comparing
  against field-measured stems. Now matches `rasterio.transform.xy`, verified
  on a whole-pixel and a subpixel top.

  Polygons were never affected: `save_segments()` goes through
  `rasterio.features.shapes`, which handles the convention itself. So points
  and polygons disagreed with each other.
- **`delineate()` under-reported the crown count by one when the scene had no
  background.** It computed `len(np.unique(crowns)) - 1`, subtracting for a
  label 0 that is absent whenever `mask_thresh` sits below the CHM's minimum.
  Two crowns were announced as one. The returned array was always correct;
  only the printed number was wrong.


## [0.3.0] — 2026-07-22

### Changed
- **`max_iters` now defaults to `None`, meaning natural termination.** It was
  200, and the bound bit: on `chm_33_2012.tif` 332 of 492 crowns stopped there
  with candidates still queued, and the crown count read **132 against 63**
  once the cap was lifted. More than a factor of two in the headline number,
  decided by a constant rather than by the canopy.

  The boundaries barely moved — 2.9 % of the partition — because a truncated
  grow blocks merges rather than misplacing pixels. That is what made it hard
  to notice: the segmentation looked right and the tree count did not.

  The cap protected nothing. Growth is bounded anyway, since every iteration
  either accepts a region or records a rejection and there are finitely many
  of both, and natural termination needed at most 484 iterations on that scene
  while running *faster* — 0.67 s against 2.34 s, because twice as many
  surviving crowns cost more in conflict arbitration than the extra merges
  cost in growing.

  `delineate_crowns()` had its own default of 200 and now passes `None` too.
  An explicit `max_iters` still works and warns when it binds; the message
  carries no seed id, so Python's default filter collapses what would
  otherwise have been 332 warnings into one.

## [0.2.0] — 2026-07-22

### Changed
- **An even smoothing or detection window is now refused.** An even window
  has no centre pixel, so `scipy.ndimage` places it half a pixel off and the
  result depends on which way the raster happens to be oriented. Measured:
  smoothing a 40x55 scene and its mirror image differed by up to 8.8 m at
  `ws=4`; `detect_tops` on `chm_150_2023.tif` found 397 tops and 400 on the
  mirror at `ws=4`, and 188 against 206 at `ws=6` — a 9 per cent difference
  in the tree count from orientation alone.

  `ws` was validated only as `>= 1`, and nothing in the output shows that it
  happened, so this raises rather than warns. `method="gaussian"` is exempt:
  `ws` only scales sigma there and the kernel stays symmetric.

  The same guard went into rcacumen, which reproduced the asymmetry faithfully.


All notable changes to this project are documented here.
The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and the project adheres to [Semantic Versioning](https://semver.org/).

## [0.1.0] — 2026-07-15

First release. pycacumen supersedes `pycrown_simplified`, which is archived.

### Added
- `CrownDelineator` — chainable pipeline: smooth → detect → merge → screen → delineate.
- `delineate_crowns()` — the whole pipeline in one call.
- Array-first API: `smooth_chm`, `detect_tops`, `merge_tops`, `screen_tops`,
  `as_pixels`, `HierarchicalRegionGrower`, `resolve_conflicts`.
- `conflict_rule` on `run_all()` — `'height'`, `'distance'` or `'similarity'`
  decide which crown wins canopy claimed by two. Ties break to the lower
  crown id, so the segmentation is reproducible whatever order the tree
  tops arrive in.
- `protect_seeds` — guarantee one crown per tree top, disabling merging.
- `retry_rejected` — reconsider regions rejected earlier in a grow.
- `n_contested` — how many regions were claimed by more than one crown.
- `read_chm()` with windowed reads, for rasters that do not fit in memory.
- Command line: `pycacumen -i chm.tif -o crowns.tif`, and `python -m pycacumen`.
- Test suite (65 tests) and CI.

### Changed
- **Split from PyCrown's structure.** Hierarchical region growing is now the
  only delineation method (see *Removed*). Smoothing, tree top detection and
  the algorithm each live in their own module, and none of them touch the disk
  — file access is confined to `pycacumen.io`.
- **Contested regions are arbitrated, not overwritten.** Previously the crown
  with the higher id silently won, so which of two merging trees survived
  depended on tree-top ordering.
- Statistics merging is attributed correctly to the pairwise formula of
  Chan, Golub & LeVeque (1983). The code called it Welford's method, which is
  a different algorithm — it adds one sample at a time, whereas a crown
  absorbs a whole region at once.
- Parallel growing shares the graph per worker instead of per task. The
  CSR graph was previously pickled once per seed: 122 MB of IPC for 1156
  trees, versus 0.4 MB now. Results are identical to sequential.
- `tree_detection` (now `detect_tops`) returns `(n, 2)` float64 rather than a
  list of tuples, matching every other tree-top function.

### Fixed
- **Growing crashed on any CHM without trees.** An empty detection produced a
  `(0,)` array instead of `(0, 2)`, so delineation raised
  `IndexError: too many indices for array`. Tiles covering clearcuts, gaps or
  water hit this.
- **Growing was unusable on in-memory arrays.** The algorithm demanded a file
  path, so `PyCrown(chm_array=...)` raised
  `TypeError: invalid path or file: None` on `hierarchical_crown_delineation`.
- Tree top correction returned early with fewer than two tops without storing
  its result, leaving the object's state inconsistent with every other path.
- Tree tops outside the raster are rejected explicitly; negative coordinates
  previously wrapped around silently.
- Dead code removed: an unused inverse-statistics helper, an unused loop
  variable in the merge formula, and a branch that could never be taken.

### Packaging
- Licence declared as an SPDX expression with `license-files` (PEP 639)
  rather than a TOML table, which setuptools deprecated and will stop
  accepting on 2027-02-18. The redundant licence classifier is gone.
  Requires setuptools >= 77 to build.
- Packages resolved with an explicit `find` directive rather than a
  hard-coded list, so `tests/` and `examples/` are never shipped and any
  future subpackage is picked up instead of being silently dropped.
  Reported against `pycrown_simplified` by Bas van Driel.

### Removed
- **Dalponte & Coomes delineation.** It is well covered by
  [PyCrown](https://github.com/manaakiwhenua/pycrown),
  [lidR](https://github.com/r-lidar/lidR) and
  [itcSegment](https://cran.r-project.org/package=itcSegment); duplicating it
  here served no one. pycacumen does one thing.
- Point cloud I/O. pycacumen takes a CHM raster and returns crowns.
