# rcacumen 0.6.1

## Changed

* **The chain warning of `merge_tops()` is written for someone who has not read the
  code.** It speaks of chains of tops and no longer of groups. Every number comes with
  its unit: chains "longer than 25 px", where "more than 5 times the merge distance"
  stood beside a distance of "5 px". It points to the spacing of the detected tops and
  to the detection window, where it named the spacing of the trees: tops are detected
  closer together than the trees stand, the more so the smaller the window. And it has
  a sentence of its own for the case when all the tops were merged into one. The text
  is the same in pycacumen 0.7.1.
  When the warning is given has not changed, and no result has.


# rcacumen 0.6.0

## Changed

* **A top goes to the nearest cell.** A top is a position in array coordinates, where
  whole numbers are the centres of the cells, and the cell under it was found by
  cutting the coordinate down: 187.67 went to cell 187, though it stands in 188. For a
  top that is the centre of a plateau, nearly half of the tops on the tiles of
  `test_data`, that put the seed one cell up or left of the centre, and now and then on
  the ground beside the crown, where a crown then grew. `as_pixels()`, `screen_tops()`,
  `grow_crowns()` and the fall-back of `merge_tops()` use the nearest cell now. A half
  goes up, and so does anything within a millionth of a pixel of one: the centre of a
  plateau two cells long is a half that comes out with a last bit either way, and that
  bit must not choose the cell. Tops and crowns change wherever a seed moves: the README's call on
  `chm_150_2014.tif` gives 215 tops and 214 crowns where 0.5.3 gave 213 and 212.
  With the seeds on their plateaus the two packages also agree more closely. All
  seventy synthetic scenes are identical with pycacumen 0.7.0, where the ten large ones
  differed on 0.032 % of pixels, and on the seven real tiles the finished crowns differ
  on 0.10 % of pixels, where they differed on 0.64 %. The crowns are steadier as well:
  a nanometre of noise on a tile moves up to 24 % of its canopy pixels to another
  crown, where it moved up to 85 % (`tools/measure_tie_sensitivity.py`).
* **A merged group whose centroid stands where no top could be is represented by its
  highest top**, and not only a group whose centroid has no data. `merge_tops()` takes
  `hmin`, the height the tops were detected above. With it and the CHM, a centroid on a
  cell that is not above `hmin`, the ground in the gap between two crowns for one, gives
  way to the highest top of the group. `delineate()` passes its own `hmin`.
  `merge_tops()` refuses a `distance` below zero and an `hmin` that is not a number.

## Added

* **`merge_tops()` warns when it merges in chains.** Grouping is transitive, and with a
  distance near the spacing of the trees a whole stand folds into a few tops: 1,264
  tops into one at 8 px on `chm_33_2012.tif`. The warning says how many groups stretch
  over more than five times the merge distance and how long the longest is. The result
  is unchanged.


# rcacumen 0.5.3

## Changed

* **Tops come in the order pycacumen lists them, and so do the crown ids.**
  `detect_tops()` listed the tops in the order its labelling met them, column
  by column, and `merge_tops()` listed the groups by an internal id; pycacumen
  reads the raster row by row and lists a group where its first member stood.
  The tops were the same and their order was not, and a crown id is the row
  number of its top: crown 7 of one package was another tree in the other.
  Both functions follow pycacumen's order now, so the same CHM gives the same
  ids in both packages wherever the crowns agree. On the seven tiles of
  `test_data` the crowns themselves are unchanged, pixel for pixel; only their
  numbers differ from 0.5.2. A crown can change where two tops tie exactly:
  a region claimed by two equally tall trees under `conflict_rule = "height"`,
  or a merged group on a hole with two equally high tops, goes to the top that
  comes first, and that is now the one it goes to in pycacumen.

## Fixed

* **`detect_tops()` and `smooth_chm()` stopped with "subscript out of bounds"
  on a raster narrower than half the window**, a single row or column with
  `ws = 5` for one. The reflected padding now folds as often as it has to, as
  scipy's does, and the results on such rasters are pycacumen's.
* An unused variable in `src/rag.c`.


# rcacumen 0.5.2

## Fixed

* **A group of merged tops whose centroid fell on a hole left no crown.**
  `merge_tops()` replaces tops closer than the radius by their centroid, and on
  a CHM with holes the centroid can fall on a cell with no data. Since 0.5.1
  such a top is dropped by `screen_tops()` instead of stopping the run, but the
  trees it stood for got no crown. `merge_tops()` takes the CHM as an optional
  `chm` argument now and represents such a group by the highest of its own
  tops; `delineate()` passes it. On the tiles of `test_data` with the README's
  parameters one tree comes back on each of `chm_150_2023.tif` (214 to 215
  tops) and `chm_33_2012.tif` (391 to 392); the other five are unchanged, and
  pycacumen 0.6.3 finds the same tops. `merge_tops(tops, distance)` without the
  CHM does what it did.


# rcacumen 0.5.1

## Fixed

* **`delineate()` stopped with "missing value where TRUE/FALSE needed" on a
  CHM with holes.** With `merge_distance` and `screen_hmin` both set, as in the
  README's quick start, a merged top could land on a cell with no data:
  `merge_tops()` puts the centroid of a group wherever it falls. `screen_tops()`
  compared that cell with `hmin`, got `NA`, and an `NA` in a logical index
  keeps the row as `(NA, NA)` instead of dropping it; `grow_crowns()` then
  stopped on it. Two of the seven real CHMs in `test_data` did this
  (`chm_150_2023.tif`, `chm_33_2012.tif`). Such a top is now dropped, which is
  what pycacumen does, because `NaN >= x` is false there. That is agreement,
  not a good answer: the tops merged into it are gone with it, and in both
  packages the trees they marked get no crown. One merged top goes this way
  on each of the two tiles.

  On the seven tiles the two packages now find the same number of tops and
  of crowns: 152, 212, 214, 386, 39, 165 and 145 crowns.

  `tools/quick_start_on_test_tiles.R` runs the README's call on every tile
  of the repository and holds the number of tops against pycacumen's, and
  the parity workflow runs it. The synthetic suite never called
  `merge_tops()` or `screen_tops()`, which is how a documented call could
  fail on documented data.

* `screen_tops()` refuses an `hmin` that is `NA` instead of returning rows
  of `NA`, and drops a top whose own coordinate is `NA`.

## Changed

* `?read_chm` said nodata is left as it is. It is not: cells the file
  declares as nodata arrive as `NA`, as terra has always read them. The page
  now says so, says that smoothing gives the rim of a hole a height, and
  says what happens to a nodata value the file does not declare.
  `?smooth_chm` still said `NA` spreads; it has been skipped since 0.4.0.

* The README's quick start ended in `plot(crowns_to_polygons(res$crowns))`,
  which fails with only rcacumen attached, because `plot` for a `SpatVector`
  is terra's. It now reads `terra::plot(...)`.

* The agreement with pycacumen on the synthetic suite is stated as measured:
  identical on its sixty small scenes, within 0.032 % of pixels on its ten
  large ones. The headline sentences said identical throughout; the table
  below them never did.


# rcacumen 0.5.0

## Changed

* **Renamed from `rHRG`.** The package name is now `rcacumen`, which is a
  breaking change for `library()` calls and for any dependency on it. The
  algorithms and the results are unchanged.


# rcacumen 0.4.2

## Changed

* `src/init.c` now calls `R_forceSymbols(dll, TRUE)` alongside
  `R_useDynamicSymbols(dll, FALSE)`, matching rgeoadaptels. The native
  routines are no longer reachable by name, so `.Call("C_rag", ...)` fails
  with *symbol name "C_rag" not in load table* instead of resolving. Both
  call sites in `R/hrg.R` already pass the registered object, so nothing
  inside the package changes; what changes is that a mistyped or renamed
  entry point fails at load rather than in the middle of someone's run, and
  the C routines stop being callable from outside rcacumen.


# rcacumen 0.4.1

## Fixed

* `man/grow_crowns.Rd` still documented `max_iters = 200L` after 0.3.0 changed
  the default to `NULL`. `R CMD check` reported it as a codoc mismatch and the
  R-CMD-check workflow, which treats warnings as failures, went red on all five
  platforms. Anyone reading `?grow_crowns` was told the cap was still there.

* The cause was structural, so that is what was fixed. Every file in `man/`
  carried a hand-written header instead of roxygen2's marker, so
  `roxygenise()` refused to touch any of them and the documentation had been
  maintained by hand — free to drift from the roxygen comments it claimed to
  come from, silently, except where `codoc` happens to look. `man/` is now
  generated, and `Roxygen: list(markdown = TRUE)` is set so the markdown in
  those comments renders instead of reaching the page as `[as_pixels()]`.

  The regenerated pages are strictly richer than the ones they replace (1813
  words against 1639); nothing was lost.


# rcacumen 0.4.0

## Fixed

* **`delineate()` now works on a real CHM file with nodata.** It did not.
  `read_chm()` returns `NA` where the file declares nodata -- terra honours the
  raster's NA flag, whatever the old docstring claimed -- `smooth_chm()` then
  propagated it and *grew* the hole, and `detect_tops()` stopped with
  "missing value where TRUE/FALSE needed", which points nowhere. The
  documented one-call entry point could not process the package's own test
  data.

  `NA` is now skipped rather than propagated: window statistics are taken over
  the cells that exist, and a window holding nothing but `NA` stays `NA`.
  `.prepare_mask()` treats `NA` as outside the mask, which is where numpy
  arrives by a different route, since `NaN > x` is `FALSE`.

  On `chm_150_2023.tif` R and pycacumen now both return **255 crowns**. Neither
  managed it before: R crashed and Python silently smoothed 3 % of the raster
  to arbitrary values, because scipy's comparison filters are undefined on
  NaN rather than NaN-aware. Matching Python would have meant reproducing that
  -- the same trap that gave this package its watershed divergence -- so both
  were defined instead. pycacumen 0.5.0 carries the other half.

  The shared cross-check gained ten scenes with nodata holes, so this is
  covered rather than assumed.

# rcacumen 0.3.0

* **`max_iters` now defaults to `NULL`, meaning natural termination.** It was
  200, and the bound bit: on `chm_33_2012.tif` 332 of 492 crowns stopped there
  with candidates still queued, and the crown count read **132 against 63**
  once the cap was lifted -- more than a factor of two in the headline number,
  decided by a constant rather than by the canopy.

  The boundaries barely moved, 2.9 % of the partition, because a truncated
  grow blocks merges rather than misplacing pixels. That is what made it hard
  to notice: the segmentation looked right and the tree count did not.

  The cap bought nothing. Growth is bounded anyway -- each iteration either
  accepts a region or records a rejection, and there are finitely many of both
  -- and natural termination needed at most 484 iterations on that scene while
  running *faster*, since twice as many surviving crowns cost more in conflict
  arbitration than the extra merges cost in growing.

  Passing an explicit `max_iters` still works and now warns when it binds.

  pycacumen 0.3.0 carries the same change.

# rcacumen 0.2.0

* **An even smoothing or detection window is now refused.** An even window has
  no centre pixel, so it sits half a pixel off and the result depends on which
  way the raster happens to be oriented. Measured before the guard: smoothing a
  40x55 scene and its mirror image differed by up to 10.5 m at `ws = 4`, and
  `detect_tops()` found 151 tops against 137 on the mirror. `method =
  "gaussian"` is exempt, because `ws` only sets sigma there and the kernel stays
  symmetric.

  pycacumen 0.2.0 carries the same guard; the asymmetry was faithfully reproduced
  here, so it had to be fixed on both sides to keep the two in step.

# rcacumen 0.1.0

First release.

* `delineate()` runs the whole pipeline; `smooth_chm()`, `detect_tops()`,
  `merge_tops()`, `screen_tops()` and `grow_crowns()` are each usable alone.
* `conflict_rule` decides which crown wins canopy claimed by two: `"height"`,
  `"distance"` or `"similarity"`. Ties resolve to the lower crown id, so the
  segmentation does not depend on the order the tree tops arrive in.
* `protect_seeds` guarantees one crown per tree top, disabling merging.
* `retry_rejected` reconsiders regions rejected earlier in a grow.
* `read_chm()`, `write_crowns()` and `crowns_to_polygons()` bridge to
  \pkg{terra}, which is only a Suggests: the algorithm works on plain matrices
  with no spatial dependency.
* Verified identical to the pycacumen Python package on 60 random scenes
  (120,977 pixels, no difference), and the watershed alone on a further 200
  scenes against `scikit-image`. See `tools/cross_validate_against_pycacumen.R`.
* On real canopy height models the agreement is close but not exact: 375 of
  147,368 masked pixels differ in the watershed (0.25 %), which the region
  growing amplifies to 5,142 of 281,602 crown pixels. `src/watershed.c`
  reimplements `skimage.segmentation.watershed` rather than calling it, and
  the two break plateau ties differently. Documented in the README rather
  than described as equality.
* `tools/measure_tie_sensitivity.py` quantifies how much the method itself
  depends on those ties. A nanometre of noise removes 3 to 30 tree tops and
  moves crowns for up to 62 % of canopy pixels — but 0 % once the tops are
  held fixed, so the instability is in tree detection, not in the growing.
