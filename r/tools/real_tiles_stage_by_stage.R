# Where rcacumen and pycacumen part on real canopy height models, stage by stage.
#
# The synthetic scenes of tools/cross_validate_against_pycacumen.R agree
# exactly or nearly so; real canopy height models do not, because they are
# full of plateaus and the two watersheds break ties differently. This puts
# the numbers on that which the README quotes, on the tiles of the
# repository, and pins the difference to the stage that makes it: every
# stage of rcacumen is fed pycacumen's own output of the stage before.
#
#   smoothing   rcacumen on the heights as read
#   detection   rcacumen on pycacumen's smoothed model
#   mask        rcacumen on pycacumen's smoothed model
#   watershed   rcacumen's on pycacumen's smoothed model, tops and mask,
#               against scikit-image's
#   pipeline    each package on its own, from the heights as read
#
# The first three have to agree and fail the run when they do not. The
# watershed and the pipeline are reported: what they differ by is the known
# plateau divergence, and its size depends on the version of scikit-image.
#
# Settings are those of tools/measure_tie_sensitivity.py: median smoothing
# over 3 px, tops above 5 m in a 5 px window, variance threshold 8, mask at
# 1 m, conflicts to the taller top; no-data is read as 0 m.
#
#   pip install "../python[io]"
#   python3 tools/generate_pycacumen_reference.py
#   Rscript tools/real_tiles_stage_by_stage.R          # from r/
#
# Copyright (C) 2026 Igor Pawelec. Licence: GPLv3.

library(rcacumen)

dir <- file.path("tools", "reference")
idx_file <- file.path(dir, "real.csv")
if (!file.exists(idx_file))
  stop("no ", idx_file, ": run tools/generate_pycacumen_reference.py first, ",
       "with rasterio installed (pip install \"../python[io]\")", call. = FALSE)
idx <- utils::read.csv(idx_file, stringsAsFactors = FALSE)
if (!nrow(idx)) stop(idx_file, " lists no tiles", call. = FALSE)

HMIN <- 5; WS <- 5L; SMOOTH_WS <- 3L; VT <- 8; MASK <- 1; RULE <- "height"

rd <- function(stem, what) {
  m <- as.matrix(utils::read.csv(
    file.path(dir, sprintf("real_%s_%s.csv", stem, what)), header = FALSE,
    na.strings = c("nan", "NA", "NaN")))
  dimnames(m) <- NULL
  m
}

# Crown ids are the row numbers of the tops a package grew from, and the two
# packages need not list their tops in one order. A pixel is compared by the
# position of the top its crown belongs to, not by the id.
seed_of <- function(out, tops) rbind(c(-1L, -1L), tops)[as.vector(out) + 1L, , drop = FALSE]

fail <- character(0)
tot <- c(mask = 0, ws = 0, px = 0, out = 0, canopy = 0)
worst <- list(tile = NA_character_, ws = 0, out = 0, ratio = -1)
sm_dev <- 0; top_dev <- 0

cat(sprintf("%-16s %5s %10s %10s %9s %10s %9s\n", "tile", "tops", "masked px",
            "watershed", "pixels", "pipeline", "of canopy"))
cat(strrep("-", 76), "\n")
for (i in seq_len(nrow(idx))) {
  stem <- idx$tile[i]
  raw <- rd(stem, "chm"); sm_py <- rd(stem, "sm"); sub <- rd(stem, "sub")
  mask_py <- rd(stem, "mask") != 0
  ws_py <- rd(stem, "ws"); out_py <- rd(stem, "out")
  if (ncol(sub) == 1L) sub <- t(sub)
  nr <- nrow(sm_py); nc <- ncol(sm_py)

  # 1. smoothing
  sm_r <- smooth_chm(raw, ws = SMOOTH_WS, method = "median")
  d_sm <- max(abs(sm_r - sm_py))
  sm_dev <- max(sm_dev, d_sm)
  if (d_sm > 0)
    fail <- c(fail, sprintf("%s: smoothed model differs by up to %.3e m", stem, d_sm))

  # 2. detection
  own <- as.matrix(detect_tops(sm_py, hmin = HMIN, ws = WS))
  if (nrow(own) != nrow(sub)) {
    fail <- c(fail, sprintf("%s: %d tops, pycacumen has %d", stem, nrow(own), nrow(sub)))
  } else {
    a <- own[order(own[, 1], own[, 2]), , drop = FALSE]
    b <- sub[order(sub[, 1], sub[, 2]), , drop = FALSE]
    top_dev <- max(top_dev, max(abs(a - b)))
    # Close is not enough at a pixel boundary: a top that floors to another
    # pixel grows another crown.
    if (any(floor(a) != floor(b)))
      fail <- c(fail, sprintf("%s: %d top(s) fall in another pixel than pycacumen's",
                              stem, sum(rowSums(floor(a) != floor(b)) > 0)))
  }

  # 3. mask
  mask_r <- rcacumen:::.prepare_mask(sm_py, MASK, 0L) != 0
  if (any(mask_r != mask_py))
    fail <- c(fail, sprintf("%s: the masks differ in %d px", stem, sum(mask_r != mask_py)))

  # 4. watershed, from the same markers on the same surface
  tops_py <- matrix(as.integer(floor(sub)), ncol = 2)
  markers <- matrix(0L, nr, nc)
  markers[tops_py] <- seq_len(nrow(tops_py))
  ws_r <- rcacumen:::.watershed(-sm_py, markers, matrix(as.integer(mask_py), nr, nc))
  d_ws <- sum(ws_r != ws_py & mask_py)
  n_mask <- sum(mask_py)

  # 5. the whole pipeline, each package on its own
  tops_r <- as_pixels(detect_tops(sm_r, hmin = HMIN, ws = WS))
  out_r <- grow_crowns(sm_r, tops_r, variance_thresh = VT, mask_thresh = MASK,
                       conflict_rule = RULE)
  moved <- rowSums(seed_of(out_r, tops_r) != seed_of(out_py, tops_py)) > 0
  d_out <- sum(moved)
  canopy <- sum(as.vector(out_r) > 0 | as.vector(out_py) > 0)

  cat(sprintf("%-16s %5d %10d %10d %9d %10d %8.2f%%\n", stem, nrow(tops_py),
              n_mask, d_ws, length(out_py), d_out, 100 * d_out / max(canopy, 1)))
  tot <- tot + c(n_mask, d_ws, length(out_py), d_out, canopy)
  ratio <- if (d_ws > 0) d_out / d_ws else -1
  if (ratio > worst$ratio)
    worst <- list(tile = stem, ws = d_ws, out = d_out, ratio = ratio)
}

big <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
cat(strrep("-", 76), "\n")
cat(sprintf("tiles                         : %d\n", nrow(idx)))
cat(sprintf("smoothing, max difference     : %.3e m\n", sm_dev))
cat(sprintf("detection, max top deviation  : %.3e px\n", top_dev))
cat(sprintf("watershed vs scikit-image     : %s of %s masked px differ (%.2f %%)\n",
            big(tot[["ws"]]), big(tot[["mask"]]), 100 * tot[["ws"]] / tot[["mask"]]))
cat(sprintf("whole pipeline vs pycacumen   : %s of %s px differ (%.2f %%), %.2f %% of the %s canopy px\n",
            big(tot[["out"]]), big(tot[["px"]]), 100 * tot[["out"]] / tot[["px"]],
            100 * tot[["out"]] / tot[["canopy"]], big(tot[["canopy"]])))
if (!is.na(worst$tile))
  cat(sprintf("largest amplification         : %s, %s watershed px become %s crown px\n",
              worst$tile, big(worst$ws), big(worst$out)))

if (top_dev >= 1e-9)
  fail <- c(fail, sprintf("tree tops deviate by %.3e px", top_dev))
if (length(fail))
  stop("rcacumen and pycacumen disagree before the watershed:\n  ",
       paste(fail, collapse = "\n  "), call. = FALSE)
cat("\nsmoothing, detection and mask agree on every tile; the difference enters\n")
cat("at the watershed and the region growing carries it further.\n")
