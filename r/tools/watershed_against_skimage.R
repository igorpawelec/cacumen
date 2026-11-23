# The watershed of rcacumen against scikit-image's, on random scenes.
#
# src/watershed.c is a reimplementation of skimage.segmentation.watershed,
# not a call to it, and everything the grower does rests on its regions. This
# holds the two against each other pixel by pixel on 200 small random
# surfaces, two in three with a mask, every second one rounded to whole
# metres so that it is plateaus from edge to edge. Each surface is flooded
# twice:
#
#   from the tops   markers on the tops the package detects, which is how
#                   the method uses the watershed;
#   from drops      the same number of markers dropped at random inside the
#                   mask, so that several of them meet on one plateau and
#                   the order of flooding decides who gets it.
#
# The README quotes what this prints. Flooded from the tops, the smooth
# surfaces have to come out identical, and the run fails when they do not.
# The rest is reported: what the two watersheds differ by there is the known
# difference in breaking plateau ties, the one
# tools/real_tiles_stage_by_stage.R measures on real canopy height models.
#
#   pip install "../python[io]"
#   python3 tools/generate_pycacumen_reference.py
#   Rscript tools/watershed_against_skimage.R          # from r/
#
# Copyright (C) 2026 Igor Pawelec. Licence: GPLv3.

library(rcacumen)

f <- file.path("tools", "reference", "watershed_scenes.csv")
if (!file.exists(f))
  stop("no ", f, ": run tools/generate_pycacumen_reference.py first", call. = FALSE)
d <- utils::read.csv(f)

kinds <- c("smooth", "rounded"); floods <- c("top", "drop")
res <- array(0, c(2, 2, 4), list(kinds, floods, c("scenes", "px", "scenes_diff", "px_diff")))
for (k in sort(unique(d$scene))) {
  s <- d[d$scene == k, ]
  nr <- max(s$row); nc <- max(s$col)
  as_matrix <- function(v) { m <- matrix(0, nr, nc); m[cbind(s$row, s$col)] <- v; m }
  height <- as_matrix(s$height); mask <- as_matrix(s$mask)
  kind <- if (all(height == round(height))) "rounded" else "smooth"
  for (fl in floods) {
    got <- rcacumen:::.watershed(-height, as_matrix(s[[paste0("marker_", fl)]]), mask)
    dk <- sum(got != as_matrix(s[[paste0("label_", fl)]]))
    res[kind, fl, ] <- res[kind, fl, ] + c(1, nr * nc, dk > 0, dk)
  }
}

big <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
cat(sprintf("%-18s %-8s %7s %10s %12s %9s\n", "markers", "surface", "scenes", "px",
            "px differing", "on scenes"))
cat(strrep("-", 70), "\n")
for (fl in floods) for (kind in kinds)
  cat(sprintf("%-18s %-8s %7d %10s %12s %9d\n",
              if (fl == "top") "on detected tops" else "dropped at random", kind,
              res[kind, fl, "scenes"], big(res[kind, fl, "px"]),
              big(res[kind, fl, "px_diff"]), res[kind, fl, "scenes_diff"]))
cat(strrep("-", 70), "\n")
for (fl in floods)
  cat(sprintf("%-18s all: %s of %s px differ (%.2f %%), on %d of %d scenes\n",
              if (fl == "top") "on detected tops" else "dropped at random",
              big(sum(res[, fl, "px_diff"])), big(sum(res[, fl, "px"])),
              100 * sum(res[, fl, "px_diff"]) / sum(res[, fl, "px"]),
              sum(res[, fl, "scenes_diff"]), sum(res[, fl, "scenes"])))

if (res["smooth", "top", "px_diff"] > 0)
  stop("flooded from the detected tops, the watershed of rcacumen differs from ",
       "scikit-image's on ", res["smooth", "top", "scenes_diff"],
       " smooth scene(s), where the two have to agree", call. = FALSE)
cat("\nflooded from the detected tops, the two watersheds agree on every smooth scene\n")
