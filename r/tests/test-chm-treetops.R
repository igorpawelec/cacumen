# Plain-R tests: any error fails R CMD check. No testthat dependency, so
# these run anywhere the package installs.
library(rcacumen)

ok <- function(cond, what) {
  if (!isTRUE(cond)) stop("FAILED: ", what, call. = FALSE)
  cat("  ok:", what, "\n")
}

cat("smoothing\n")
m <- matrix(1, 5, 5); m[3, 3] <- 99
ok(smooth_chm(m, 3, "median")[3, 3] == 1, "median kills a spike")
ok(identical(dim(smooth_chm(m, 3)), dim(m)), "shape preserved")
ok(all(smooth_chm(m, 1) == m), "ws=1 is a no-op")
for (meth in c("median", "mean", "gaussian", "maximum"))
  ok(all(is.finite(smooth_chm(m, 3, meth))), paste(meth, "runs"))
ok(inherits(try(smooth_chm(m, 0), silent = TRUE), "try-error"), "ws=0 rejected")
ok(inherits(try(smooth_chm(m, 3, "nope"), silent = TRUE), "try-error"),
   "bad method rejected")

cat("tree tops\n")
chm <- matrix(0, 30, 30); chm[10, 10] <- 25; chm[22, 20] <- 18
sm <- smooth_chm(chm, 3, "maximum")
t <- detect_tops(sm, hmin = 5, ws = 5)
ok(nrow(t) == 2, "two peaks found")
ok(ncol(t) == 2, "tops are (row, col)")
ok(is.matrix(t) && is.numeric(t), "tops are a numeric matrix")

e <- detect_tops(matrix(0, 10, 10), hmin = 5)
ok(identical(dim(e), c(0L, 2L)), "empty detection is 0x2, not 0")
invisible(e[, 1])   # must not error

ok(nrow(detect_tops(sm, hmin = 100, ws = 5)) == 0, "hmin filters everything")

cat("merging\n")
tp <- rbind(c(10, 10), c(10, 12), c(50, 50))
ok(nrow(merge_tops(tp, 5)) == 2, "close pair merges")
ok(nrow(merge_tops(tp, 0.5)) == 3, "distant tops kept")
ok(all(merge_tops(rbind(c(10, 10), c(10, 20)), 20)[1, ] == c(10, 15)),
   "merged top is the centroid")
ok(identical(dim(merge_tops(matrix(numeric(0), 0, 2))), c(0L, 2L)), "empty merge")

cat("screening\n")
ok(nrow(screen_tops(sm, t, 100)) == 0, "screen drops short tops")
ok(nrow(screen_tops(sm, t, 0)) == 2, "screen keeps tall tops")
ok(nrow(screen_tops(sm, rbind(c(-5, -5), c(500, 500)), 0)) == 0,
   "out-of-raster tops dropped, not wrapped")

# A top on a cell with no data used to survive as a row of NA, on which
# grow_crowns() stopped with "missing value where TRUE/FALSE needed". It
# happened on two of the seven real CHMs in test_data with the README's own
# parameters: merge_tops() puts the centroid of a group wherever it falls.
holed <- matrix(20, 10, 10); holed[1:3, ] <- NA
kept <- screen_tops(holed, rbind(c(2, 5), c(7, 5)), hmin = 10)
ok(nrow(kept) == 1 && !anyNA(kept) && all(kept[1, ] == c(7, 5)),
   "a top on an NA cell is dropped, not kept as an NA row")
# The same NA rows by two other doors: an hmin that is NA, and a top whose
# own coordinate is NA.
ok(inherits(try(screen_tops(holed, rbind(c(7, 5)), NA), silent = TRUE), "try-error"),
   "hmin = NA is refused, not turned into NA rows")
kept <- screen_tops(holed, rbind(c(NA, 1), c(7, 5)), hmin = 10)
ok(nrow(kept) == 1 && all(kept[1, ] == c(7, 5)),
   "a top with an NA coordinate is dropped")

# The same through delineate(): two crowns either side of a hole, close
# enough to merge, so their centroid lands inside the hole. Until 0.5.2 the
# merged top stayed there and was screened out, and both trees were left
# without a crown; now the higher of the two stands in for the group.
yy <- row(matrix(0, 30, 30)); xx <- col(matrix(0, 30, 30))
gap <- pmax(26 * exp(-((yy - 15)^2 + (xx - 10)^2) / 8),
            26 * exp(-((yy - 15)^2 + (xx - 16)^2) / 8))
gap[13:17, 12:14] <- NA
res <- delineate(gap, smooth_ws = 1, hmin = 5, detect_ws = 5,
                 merge_distance = 8, screen_hmin = 10, quiet = TRUE)
ok(!anyNA(res$tops), "delineate() keeps no NA top when a merged top falls in a hole")
ok(nrow(res$tops) == 1 &&
   !is.na(gap[floor(res$tops[1, 1]), floor(res$tops[1, 2])]) &&
   sum(res$crowns > 0) > 0,
   "the group whose centroid falls in the hole keeps a top and a crown")

cat("as_pixels\n")
ok(all(as_pixels(rbind(c(3.7, 4.2))) == c(3L, 4L)), "floors to integers")
ok(is.integer(as_pixels(rbind(c(3.7, 4.2)))), "returns integers")


# ── even windows are refused ──────────────────────────────────────────
#
# An even window has no centre pixel, so it sits half a pixel off and the
# result depends on which way the raster is oriented. Measured before the
# guard: smoothing a 40x55 scene and its mirror differed by up to 10.5 m at
# ws=4, and detect_tops found 151 tops against 137 on the mirror.
cat("even windows\n")

set.seed(4)
sc <- matrix(runif(40 * 55) * 20 + 5, 40, 55)
bad <- function(expr, what) ok(inherits(try(expr, silent = TRUE), "try-error"), what)

for (ws in c(2L, 4L, 6L)) {
  for (meth in c("median", "mean", "maximum"))
    bad(smooth_chm(sc, ws = ws, method = meth),
        sprintf("smooth_chm rejects ws=%d with method='%s'", ws, meth))
  bad(detect_tops(sc, hmin = 5, ws = ws),
      sprintf("detect_tops rejects ws=%d", ws))
}

# gaussian is exempt: ws only sets sigma, and the kernel stays symmetric.
g <- smooth_chm(sc, ws = 4L, method = "gaussian")
gm <- smooth_chm(sc[, ncol(sc):1], ws = 4L, method = "gaussian")[, ncol(sc):1]
ok(isTRUE(all.equal(g, gm)), "gaussian accepts an even ws and stays symmetric")

for (meth in c("median", "mean", "maximum", "gaussian")) {
  a <- smooth_chm(sc, ws = 5L, method = meth)
  b <- smooth_chm(sc[, ncol(sc):1], ws = 5L, method = meth)[, ncol(sc):1]
  ok(isTRUE(all.equal(a, b)),
     sprintf("odd ws is orientation-independent for '%s'", meth))
}
ok(nrow(detect_tops(sc, hmin = 5, ws = 5L)) ==
   nrow(detect_tops(sc[, ncol(sc):1], hmin = 5, ws = 5L)),
   "detection with an odd ws is orientation-independent")

cat("a merged top on a hole\n")

# Two crowns six cells apart with a hole between them, wide enough to
# outlast the smoothing, which fills the cells on a hole's rim. The
# centroid of their two tops falls on the hole.
hole_scene <- function(n = 41L) {
  rr <- matrix(seq_len(n), n, n); cc <- t(rr)
  m <- matrix(0, n, n)
  for (p in list(c(21, 18, 22), c(21, 24, 26)))
    m <- pmax(m, p[3] * exp(-((rr - p[1])^2 + (cc - p[2])^2) / (2 * 2.5^2)))
  m[19:23, 20:22] <- NA
  m
}
hs <- hole_scene()
two <- rbind(c(21.5, 18.5), c(21.5, 24.5))
ok(is.na(hs[21, 21]), "the centroid of the two tops falls on no data")
ok(isTRUE(all.equal(merge_tops(two, 7), rbind(c(21.5, 21.5)))),
   "without the CHM the centroid stands wherever it falls")
ok(isTRUE(all.equal(merge_tops(two, 7, chm = hs), rbind(c(21.5, 24.5)))),
   "with the CHM the highest top stands in for a centroid on no data")
ok(nrow(screen_tops(hs, merge_tops(two, 7, chm = hs), 10)) == 1L &&
   nrow(screen_tops(hs, merge_tops(two, 7), 10)) == 0L,
   "and survives the screening that drops the centroid")
tie <- hs; tie[21, 18] <- tie[21, 24] <- 30
ok(isTRUE(all.equal(merge_tops(two, 7, chm = tie), rbind(c(21.5, 18.5)))) &&
   isTRUE(all.equal(merge_tops(two[2:1, ], 7, chm = tie), rbind(c(21.5, 24.5)))),
   "equally high tops give the first of them")
filled <- hs; filled[19:23, 21] <- 5
ok(isTRUE(all.equal(merge_tops(two, 7, chm = filled), rbind(c(21.5, 21.5)))),
   "a centroid on data is left alone")
five <- rbind(c(6.5, 6.5), c(21.5, 18.5), c(7.5, 6.5), c(21.5, 24.5), c(36.5, 36.5))
plain <- merge_tops(five, 7); aware <- merge_tops(five, 7, chm = hs)
moved <- rowSums(plain != aware) > 0
ok(nrow(plain) == 3L && nrow(aware) == 3L && sum(moved) == 1L &&
   isTRUE(all.equal(aware[moved, ], c(21.5, 24.5))),
   "groups that are not affected keep their place and order")
ok(isTRUE(all.equal(merge_tops(rbind(c(10.5, 3.5), c(13.5, 3.5)), 5,
                               chm = matrix(20, 10, 10)), rbind(c(10.5, 3.5)))),
   "a centroid outside the raster counts as no data")
res <- delineate(hs, smooth_ws = 3, hmin = 5, detect_ws = 3, merge_distance = 7,
                 screen_hmin = 10, variance_thresh = 50, mask_thresh = 1, quiet = TRUE)
ok(nrow(res$tops) == 1L && res$tops[1, 2] > 22 && sum(res$crowns > 0) > 20,
   "delineate() keeps a crown for the group, grown from the taller tree")

cat("rasters narrower than the window\n")

# Until 0.5.3 these stopped with "subscript out of bounds": the reflected
# padding asked for more cells than a single row has. The expected values
# are what pycacumen returns for the same numbers (its coordinates plus 1).
thin <- matrix(c(1, 5, 2, 7, 3, 9, 4, 6, 2, 8, 1, 3, 2), nrow = 1)
ok(isTRUE(all.equal(detect_tops(thin, hmin = 0.5, ws = 5), rbind(c(1, 6), c(1, 10)))),
   "detect_tops() on a single row, window 5")
ok(isTRUE(all.equal(detect_tops(t(thin), hmin = 0.5, ws = 5), rbind(c(6, 1), c(10, 1)))),
   "detect_tops() on a single column, window 5")
ok(isTRUE(all.equal(as.vector(smooth_chm(thin, ws = 5, method = "median")),
                    c(2, 2, 3, 5, 4, 6, 4, 6, 4, 3, 2, 2, 2))),
   "median smoothing of a single row, window 5")
two <- rbind(c(3, 8, 1, 9, 2, 7, 4), c(6, 2, 9, 1, 5, 3, 8))
ok(isTRUE(all.equal(detect_tops(two, hmin = 0.5, ws = 5),
                    rbind(c(1, 4), c(2, 3), c(2, 7)))),
   "detect_tops() on two rows, window 5")
ok(isTRUE(all.equal(smooth_chm(two, ws = 5, method = "mean"),
                    rbind(c(4.84, 4.8, 4.6, 4.56, 4.96, 5.08, 5.16),
                          c(4.76, 4.8, 4.6, 4.84, 4.84, 5.12, 5.04)))),
   "mean smoothing of two rows, window 5")

cat("order of the tops\n")

# pycacumen holds itself to the same two orders (python/tests/test_treetops.py),
# which is what makes a crown id mean the same tree in both packages.
# Two peaks: read row by row the one in row 5 comes first, read column by
# column the one in column 5 would.
oo <- matrix(0, 30, 30); oo[5, 20] <- 20; oo[10, 5] <- 25
ok(isTRUE(all.equal(detect_tops(oo, hmin = 5, ws = 3), rbind(c(5, 20), c(10, 5)))),
   "tops come in the order of a row-by-row reading")
og <- rbind(c(10, 10), c(30, 30), c(11, 10), c(50, 50), c(31, 30))
ok(isTRUE(all.equal(merge_tops(og, 2), rbind(c(10.5, 10), c(30.5, 30), c(50, 50)))),
   "merged groups come in the order of their first member")

cat("\nall chm/treetops tests passed\n")
