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
   !is.na(gap[as_pixels(res$tops)]) &&
   sum(res$crowns > 0) > 0,
   "the group whose centroid falls in the hole keeps a top and a crown")

cat("as_pixels\n")
# Whole numbers are the centres of the cells. Until 0.6.0 the coordinate was
# cut down, and (3.7, 4.2) went to cell (3, 4).
ok(all(as_pixels(rbind(c(3.7, 4.2))) == c(4L, 4L)), "a top goes to the nearest cell")
ok(all(as_pixels(rbind(c(187.67, 49.67), c(1.49, 1), c(3, 8))) ==
       rbind(c(188L, 50L), c(1L, 1L), c(3L, 8L))), "nearest cell, more cases")
ok(all(as_pixels(rbind(c(1.5, 2.5), c(3.5, 4.5))) == rbind(c(2L, 3L), c(4L, 5L))),
   "a half goes up")
# A plateau two cells long has its centre on the boundary between them. The
# weighted mean that finds it gives 41.5 or 41.49999999999999 according to the
# order of the sums, and without a margin the two packages seeded different
# cells on a double-precision CHM. The same cases as in pycacumen's tests.
for (case in list(c(12.3, 41), c(7.7, 7), c(0.1 + 0.2, 11))) {
  pm <- matrix(0, 60, 20); pm[case[2] + 0:1, 5] <- case[1]
  pt <- detect_tops(pm, hmin = 0.05, ws = 3)
  ok(nrow(pt) == 1L && abs(pt[1, 1] - (case[2] + 0.5)) < 1e-9 &&
     all(as_pixels(pt) == c(case[2] + 1, 5)),
     sprintf("a plateau of two cells at %.1f m is seeded on its second cell", case[1]))
}
ok(all(as_pixels(rbind(c(41.5 - 1e-9, 3), c(41.5 - 1e-3, 3))) == rbind(c(42L, 3L), c(41L, 3L))),
   "the margin is a millionth of a pixel, no more")
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

# Two crowns on one row, the right-hand one the taller. With a hole, a block
# of no data between them, wide enough to outlast the smoothing, which fills
# the cells on a hole's rim. The same scenes as python/tests/test_merge_on_hole.py,
# one row and one column further on, because R counts from 1.
two_crowns <- function(n = 41L, left = 18, right = 24, hole = TRUE) {
  rr <- matrix(seq_len(n), n, n); cc <- t(rr)
  m <- matrix(0, n, n)
  for (p in list(c(21, left, 22), c(21, right, 26)))
    m <- pmax(m, p[3] * exp(-((rr - p[1])^2 + (cc - p[2])^2) / (2 * 2.5^2)))
  if (hole) m[19:23, 20:22] <- NA
  m
}
hs <- two_crowns()
two <- rbind(c(21, 18), c(21, 24))           # whole numbers: the centres of two cells
ok(is.na(hs[21, 21]), "the centroid of the two tops falls on no data")
ok(isTRUE(all.equal(merge_tops(two, 7), rbind(c(21, 21)))),
   "without the CHM the centroid stands wherever it falls")
ok(isTRUE(all.equal(merge_tops(two, 7, chm = hs), rbind(c(21, 24)))),
   "with the CHM the highest top stands in for a centroid on no data")
ok(nrow(screen_tops(hs, merge_tops(two, 7, chm = hs), 10)) == 1L &&
   nrow(screen_tops(hs, merge_tops(two, 7), 10)) == 0L,
   "and survives the screening that drops the centroid")
tie <- hs; tie[21, 18] <- tie[21, 24] <- 30
ok(isTRUE(all.equal(merge_tops(two, 7, chm = tie), rbind(c(21, 18)))) &&
   isTRUE(all.equal(merge_tops(two[2:1, ], 7, chm = tie), rbind(c(21, 24)))),
   "equally high tops give the first of them")
filled <- hs; filled[19:23, 21] <- 5
ok(isTRUE(all.equal(merge_tops(two, 7, chm = filled), rbind(c(21, 21)))),
   "a centroid on data is left alone")
five <- rbind(c(6, 6), c(21, 18), c(7, 6), c(21, 24), c(36, 36))
plain <- merge_tops(five, 7); aware <- merge_tops(five, 7, chm = hs)
moved <- rowSums(plain != aware) > 0
ok(nrow(plain) == 3L && nrow(aware) == 3L && sum(moved) == 1L &&
   isTRUE(all.equal(aware[moved, ], c(21, 24))),
   "groups that are not affected keep their place and order")
ok(isTRUE(all.equal(merge_tops(rbind(c(10, 3), c(13, 3)), 5,
                               chm = matrix(20, 10, 10)), rbind(c(10, 3)))),
   "a centroid outside the raster counts as no data")
res <- delineate(hs, smooth_ws = 3, hmin = 5, detect_ws = 3, merge_distance = 7,
                 screen_hmin = 10, variance_thresh = 50, mask_thresh = 1, quiet = TRUE)
ok(nrow(res$tops) == 1L && res$tops[1, 2] > 22 && sum(res$crowns > 0) > 20,
   "delineate() keeps a crown for the group, grown from the taller tree")

cat("a merged top in a gap\n")

# Two crowns far enough apart to leave ground between them: the centroid of
# their tops stands on that ground, where no top could have been detected.
gs <- two_crowns(left = 13, right = 29, hole = FALSE)
far2 <- rbind(c(21, 13), c(21, 29))
ok(gs[21, 21] < 0.5, "the gap is ground")
ok(isTRUE(all.equal(merge_tops(far2, 17, chm = gs), rbind(c(21, 21)))),
   "without hmin the centroid stays on the ground")
ok(isTRUE(all.equal(merge_tops(far2, 17, chm = gs, hmin = 5), rbind(c(21, 29)))),
   "with hmin the highest top stands in")
saddle <- two_crowns(left = 19, right = 23, hole = FALSE)
ok(saddle[21, 21] > 5 &&
   isTRUE(all.equal(merge_tops(rbind(c(21, 19), c(21, 23)), 5, chm = saddle, hmin = 5),
                    rbind(c(21, 21)))),
   "a centroid above hmin is left alone")
edge <- gs; edge[21, 21] <- 5
at <- merge_tops(far2, 17, chm = edge, hmin = 5)
edge[21, 21] <- 5.5
above <- merge_tops(far2, 17, chm = edge, hmin = 5)
ok(isTRUE(all.equal(at, rbind(c(21, 29)))) && isTRUE(all.equal(above, rbind(c(21, 21)))),
   "hmin is the threshold of detection itself: at it is not above it")
res <- delineate(gs, smooth_ws = 3, hmin = 5, detect_ws = 3, merge_distance = 17,
                 variance_thresh = 50, mask_thresh = 1, quiet = TRUE)
ok(nrow(res$tops) == 1L && res$tops[1, 2] > 21,
   "delineate() passes the threshold it detected with")

bad <- function(expr) inherits(try(expr, silent = TRUE), "try-error")
ok(bad(merge_tops(rbind(c(1, 1), c(1, 3)), -2)) && bad(merge_tops(rbind(c(1, 1), c(1, 3)), NA)) &&
   bad(merge_tops(rbind(c(1, 1), c(1, 3)), Inf)) && nrow(merge_tops(rbind(c(1, 1), c(1, 3)), 0)) == 2L,
   "a distance that means nothing is refused; zero merges nothing")
ok(bad(merge_tops(rbind(c(1, 1), c(1, 3)), 3, chm = matrix(1, 5, 5), hmin = NA)),
   "an hmin that is not a number is refused")

cat("chains\n")

caught <- function(expr) {
  msg <- NULL
  val <- withCallingHandlers(expr, warning = function(w) {
    msg <<- c(msg, conditionMessage(w)); invokeRestart("muffleWarning") })
  list(value = val, warnings = msg)
}
chain <- cbind(10, 4 * (0:11))               # each 4 px from the next
got <- caught(merge_tops(chain, 5))
ok(length(got$warnings) == 1L &&
   grepl("stretch over more than 5 times the merge distance", got$warnings) &&
   nrow(got$value) == 1L,
   "a merge that runs in a chain says so, and the result is what it was")
three <- rbind(chain, cbind(60, 4 * (0:7)), c(90, 10), c(90, 13))
got <- caught(merge_tops(three, 5))
ok(grepl("2 of 3 groups", got$warnings) && grepl("joins 12 tops over 44.0 px", got$warnings),
   "the warning counts the groups and names the longest")
got <- caught(merge_tops(rbind(c(10, 10), c(12, 11), c(11, 13), c(40, 40), c(42, 41)), 5))
ok(is.null(got$warnings) && nrow(got$value) == 2L, "compact groups are merged in silence")

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
