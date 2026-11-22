# The README quick start, run on every canopy height model in the repository.
#
# The quick start was written against one tile and stopped with an error on
# two of the other six: a merged top landed on a cell with no data and
# screen_tops() kept it as a row of NA. Nothing ran the documented call on
# the documented data, so nothing noticed. This does.
#
# When tools/generate_pycacumen_reference.py has been run first (with
# rasterio installed), the number of tops is also held against pycacumen's:
# the two detect, merge and screen identically, so one top more or fewer is
# a defect in one of them. Crown counts are printed side by side and
# reported when they differ, but do not fail the run: on real canopy height
# models the two watersheds break plateau ties differently, and that can
# move a merge.
#
# It lives in tools/ and runs from the parity workflow rather than from
# tests/, because test_data is .Rbuildignore'd and the tiles of the Python
# package sit outside this folder: under R CMD check neither is there.
#
#   Rscript tools/quick_start_on_test_tiles.R      # from r/
#
# Copyright (C) 2026 Igor Pawelec. Licence: GPLv3.

library(rcacumen)

dirs <- c("test_data", file.path("..", "python", "test_data"))
files <- unlist(lapply(dirs[dir.exists(dirs)], list.files,
                       pattern = "\\.tif$", full.names = TRUE))
files <- files[!duplicated(basename(files))]
files <- files[order(basename(files))]
if (!length(files))
  stop("no test tiles found in ", paste(dirs, collapse = " or "),
       "; run this from r/ in a checkout of the repository", call. = FALSE)

ref_file <- file.path("tools", "reference", "tiles.csv")
ref <- if (file.exists(ref_file)) utils::read.csv(ref_file) else NULL
if (is.null(ref))
  cat("no", ref_file, "- rcacumen is run on its own, not compared with pycacumen\n")

fail <- character(0); note <- character(0)
for (f in files) {
  tile <- basename(f)
  res <- tryCatch(
    delineate(f, smooth_ws = 3, hmin = 7, detect_ws = 5, merge_distance = 5,
              screen_hmin = 10, variance_thresh = 2, quiet = TRUE),
    error = function(e) e)
  if (inherits(res, "error")) {
    fail <- c(fail, sprintf("%s: %s", tile, conditionMessage(res)))
    next
  }
  ids <- unique(as.vector(terra::values(res$crowns)))
  n <- sum(!is.na(ids) & ids > 0)
  line <- sprintf("%-18s %4d tops, %4d crowns, %d contested",
                  tile, nrow(res$tops), n, res$n_contested)
  if (anyNA(res$tops))
    fail <- c(fail, sprintf("%s: a tree top is NA", tile))
  if (n == 0L)
    fail <- c(fail, sprintf("%s: no crowns", tile))
  if (!is.null(ref)) {
    py <- ref[ref$tile == tile, ]
    if (nrow(py) != 1L) {
      fail <- c(fail, sprintf("%s: not in %s", tile, ref_file))
    } else {
      line <- sprintf("%s | pycacumen %4d tops, %4d crowns", line, py$tops, py$crowns)
      if (py$tops != nrow(res$tops))
        fail <- c(fail, sprintf("%s: %d tops, pycacumen has %d",
                                tile, nrow(res$tops), py$tops))
      if (py$crowns != n)
        note <- c(note, sprintf("%s: %d crowns, pycacumen has %d", tile, n, py$crowns))
    }
  }
  cat(line, "\n")
}

if (length(note))
  cat("\ncrown counts differ (reported, not failed):\n  ",
      paste(note, collapse = "\n  "), "\n", sep = "")
if (length(fail))
  stop("the README quick start on the repository's own tiles:\n  ",
       paste(fail, collapse = "\n  "), call. = FALSE)
cat(sprintf("\nthe README quick start runs on all %d tiles%s\n", length(files),
            if (is.null(ref)) "" else ", with pycacumen's tops on each"))
