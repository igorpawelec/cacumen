#' Detect tree tops as local maxima of a CHM
#'
#' A pixel is a candidate when it equals the maximum of its `ws` x `ws`
#' neighbourhood and exceeds `hmin`. Adjacent candidates, which a flat apex
#' produces, are grouped into connected components and each component is
#' reduced to its height-weighted centre of mass, giving subpixel positions.
#'
#' @param chm Numeric matrix, ideally already smoothed by [smooth_chm()].
#'   On a raw CHM this over-detects badly.
#' @param hmin Numeric. Minimum height for a pixel to be considered, in the
#'   units of the CHM. Sets the floor between canopy and understorey.
#' @param ws Integer. Neighbourhood size, which acts as the minimum spacing
#'   between tops. Over-detection is the safer error: [delineate()] merges
#'   surplus tops back together, but a top never detected is a tree lost.
#'
#' @return A two-column numeric matrix of (row, col) positions, 1-based, in the
#'   order of the first cell of each maximum when the raster is read row by
#'   row, which is the order pycacumen lists them in.
#'   Zero detections give a 0 x 2 matrix, not `NULL`, so callers can index
#'   `[, 1]` unconditionally.
#'
#' @references
#' Popescu, S.C., Wynne, R.H. (2004). Seeing the trees in the forest.
#' *Photogrammetric Engineering & Remote Sensing* 70(5), 589-604.
#'
#' @examples
#' chm <- matrix(0, 20, 20); chm[10, 10] <- 25
#' detect_tops(smooth_chm(chm, 3, "maximum"), hmin = 5, ws = 3)
#' @export
detect_tops <- function(chm, hmin = 2, ws = 3L) {
  chm <- as.matrix(chm)
  if (length(dim(chm)) != 2L) stop("chm must be a 2-D matrix", call. = FALSE)
  ws <- as.integer(ws)
  if (is.na(ws) || ws < 1L) stop("ws must be >= 1", call. = FALSE)
  # An even window has no centre pixel, so the local-maximum test is taken
  # over a window sitting half a pixel off, and the tree count then depends
  # on which way the raster happens to be oriented. Measured on
  # chm_150_2023.tif: ws=4 found 397 tops and 400 on the mirror image, ws=6
  # found 188 and 206. Refused rather than warned about, because nothing in
  # the output shows it happened.
  if (ws %% 2L == 0L)
    stop("ws must be odd, got ", ws, ". An even window has no centre pixel, ",
         "so the local-maximum test sits half a pixel off and the tree count ",
         "depends on the raster's orientation. Use ", ws - 1L, " or ", ws + 1L,
         ".", call. = FALSE)

  # -Inf, not NA, for the cells with no data: it loses every comparison and
  # drops out of the maximum cleanly. A NA pixel is never detected itself,
  # since NA == anything and NA > hmin are both NA, which `&` treats as
  # false once the comparison is made explicit below.
  filled <- chm; filled[is.na(filled)] <- -Inf
  loc_max <- .focal(filled, ws, function(v) if (all(!is.finite(v))) -Inf else
                      max(v, na.rm = TRUE))
  cand <- !is.na(chm) & (chm == loc_max) & (chm > hmin)
  if (!any(cand)) return(.empty_tops())

  lab <- .label_components(cand)
  nlab <- max(lab)
  idx <- which(lab > 0L)
  g <- lab[idx]
  w <- chm[idx]
  rr <- ((idx - 1L) %% nrow(chm)) + 1L
  cc <- ((idx - 1L) %/% nrow(chm)) + 1L
  sw <- as.vector(rowsum(w, g, reorder = TRUE))
  sr <- as.vector(rowsum(w * rr, g, reorder = TRUE))
  sc <- as.vector(rowsum(w * cc, g, reorder = TRUE))
  # centre of mass; a zero-weight component falls back to its mean position
  zero <- sw == 0
  if (any(zero)) {
    sr[zero] <- as.vector(rowsum(rr, g, reorder = TRUE))[zero]
    sc[zero] <- as.vector(rowsum(cc, g, reorder = TRUE))[zero]
    sw[zero] <- as.vector(tabulate(g, nbins = nlab))[zero]
  }
  # In the order pycacumen lists them: by the first cell of each maximum,
  # reading the raster row by row. The labelling above numbers them as it
  # meets them, column by column, so without this the two packages gave the
  # same tops in different orders, and the crowns grown from them the same
  # shapes under different ids.
  first <- as.vector(tapply((rr - 1L) * ncol(chm) + cc, g, min))
  cbind(sr / sw, sc / sw)[order(first), , drop = FALSE]
}

.empty_tops <- function() matrix(numeric(0), nrow = 0, ncol = 2)

# The cell each top stands on. A top is a position in array coordinates, where
# whole numbers are the centres of the cells, so the cell under it is the
# nearest one: 187.67 stands in cell 188. Until 0.6.0 the coordinate was cut
# down, to 187, which put the seed of a plateau one cell up or left of its
# centre for nearly half of the tops, and now and then on the ground beside
# the crown.
#
# A half goes up, the same way in pycacumen, and so does anything within a
# millionth of a pixel of a half. The centre of a plateau two cells long is a
# half exactly, but it is computed as a weighted mean, and comes out as 41.5
# or as 41.49999999999999 according to the order of the sums: without the
# margin the last bit chose the cell, and chose differently in the two
# packages.
.HALF_MARGIN <- 1e-6
.cells <- function(tops) floor(tops + (0.5 + .HALF_MARGIN))

# A merged group that stretches over more than this many merge distances was
# joined as a chain, tree to tree, and merge_tops() says so.
.CHAIN_FACTOR <- 5

.as_tops <- function(x) {
  if (is.null(x) || length(x) == 0L) return(.empty_tops())
  m <- matrix(as.double(as.matrix(x)), ncol = 2)
  m
}

# 4-connected labelling, matching scipy.ndimage.label; iterative so deep
# components cannot blow the R stack.
.label_components <- function(mask) {
  nr <- nrow(mask); nc <- ncol(mask)
  lab <- matrix(0L, nr, nc)
  cur <- 0L
  idx <- which(mask)
  for (start in idx) {
    if (lab[start] != 0L) next
    cur <- cur + 1L
    stack <- start
    lab[start] <- cur
    while (length(stack)) {
      p <- stack[length(stack)]; stack <- stack[-length(stack)]
      r <- ((p - 1L) %% nr) + 1L; cc <- ((p - 1L) %/% nr) + 1L
      for (d in list(c(-1L, 0L), c(1L, 0L), c(0L, -1L), c(0L, 1L))) {
        r2 <- r + d[1]; c2 <- cc + d[2]
        if (r2 < 1L || r2 > nr || c2 < 1L || c2 > nc) next
        q <- (c2 - 1L) * nr + r2
        if (mask[q] && lab[q] == 0L) { lab[q] <- cur; stack <- c(stack, q) }
      }
    }
  }
  lab
}

#' Merge tree tops that lie close together
#'
#' Tops closer than `distance` are almost always one tree seen through
#' several local maxima. Neighbours within the radius are grouped
#' transitively and each group is replaced by its centroid.
#'
#' Grouping is transitive: a chain of tops each within `distance` of the
#' next collapses into one even if the endpoints are far apart. Keep the
#' threshold well below the typical crown diameter, and below the spacing
#' of the detected tops, which can stand as close as the detection window.
#'
#' A centroid is a point between the tops. On a CHM with holes it can fall on
#' a cell with no data, and between two crowns it can fall in the gap, on the
#' ground. Nothing grows from there, or a crown grows on the ground: the trees
#' the group stood for are lost. Given the CHM, such a group is represented by
#' the highest of its own tops instead.
#'
#' A warning is given when a chain longer than five times the merge distance
#' was merged into one top: the tops were joined tree to tree, and the
#' distance is too large for the stand. It says how many of the tops left
#' after merging stand for such a chain and how long the longest is, or that
#' all the tops were merged into one.
#'
#' @param tops Two-column matrix of (row, col) positions.
#' @param distance Numeric. Merge radius in pixels.
#' @param chm Optional numeric matrix, the CHM the tops were detected on.
#'   With it, a group whose centroid falls on a cell with no data (`NA`) or
#'   outside the raster is replaced by the highest of its tops (the first of
#'   them when several are equally high). Without it every group is replaced
#'   by its centroid, wherever that falls.
#' @param hmin Optional number, the height the tops were detected above. With
#'   it and the CHM, a group whose centroid stands on a cell that is not above
#'   `hmin`, where no top could have been detected, is replaced by its highest
#'   top as well.
#'
#' @return A two-column matrix with at most as many rows as `tops`. The
#'   groups come in the order of their first member in `tops`, as pycacumen
#'   returns them.
#' @export
merge_tops <- function(tops, distance = 5, chm = NULL, hmin = NULL) {
  if (!is.numeric(distance) || length(distance) != 1L || !is.finite(distance) || distance < 0)
    stop("distance must be a single number, zero or more", call. = FALSE)
  if (!is.null(hmin) && (!is.numeric(hmin) || length(hmin) != 1L || !is.finite(hmin)))
    stop("hmin must be a single number", call. = FALSE)
  tops <- .as_tops(tops)
  n <- nrow(tops)
  if (n < 2L) return(tops)

  parent <- seq_len(n)
  find <- function(x) { while (parent[x] != x) { parent[x] <<- parent[parent[x]]; x <- parent[x] }; x }
  d2 <- distance^2
  for (i in seq_len(n - 1L)) {
    j <- (i + 1L):n
    dd <- (tops[j, 1] - tops[i, 1])^2 + (tops[j, 2] - tops[i, 2])^2
    for (k in j[dd <= d2]) {
      ri <- find(i); rk <- find(k)
      if (ri != rk) parent[ri] <- rk
    }
  }
  g <- vapply(seq_len(n), find, integer(1))
  # Number the groups by their first member, which is the order pycacumen
  # returns them in; the root ids above depend on how the unions happened
  # to run.
  g <- match(g, unique(g))
  out <- cbind(as.vector(rowsum(tops[, 1], g, reorder = TRUE)) / as.vector(table(g)),
               as.vector(rowsum(tops[, 2], g, reorder = TRUE)) / as.vector(table(g)))
  # Grouping is transitive, so tops a merge distance apart, one after another,
  # become one group however long the row. With a distance near the spacing of
  # the detected tops a whole stand folds into a few tops, and nothing in the
  # result says that it happened. The text is read by people who have not seen
  # this code: it names no groups, gives every number with its unit, and says
  # it differently when nothing is left but one top. pycacumen gives the same
  # text.
  extent <- function(v) as.vector(tapply(v, g, function(x) max(x) - min(x)))
  spans <- sqrt(extent(tops[, 1])^2 + extent(tops[, 2])^2)
  far <- spans > .CHAIN_FACTOR * distance
  if (any(far)) {
    worst <- which.max(spans)
    what <- if (length(spans) == 1L) {
      sprintf("merge_tops: all %d tops were merged into one, a chain %.1f px long. ",
              n, spans[worst])
    } else {
      sprintf(paste0(
        "merge_tops: chains of tops longer than %g px were merged into single tops ",
        "(%d of the %d tops left after merging). The longest chain joined %d tops ",
        "over %.1f px. "),
        .CHAIN_FACTOR * distance, sum(far), length(spans),
        sum(g == worst), spans[worst])
    }
    warning(what, sprintf(paste0(
      "Tops are joined one after another, each within the merge distance (%g px) of ",
      "the next, so a distance near the spacing of the detected tops joins them ",
      "across the stand. Detected tops can stand as close as the detection window, ",
      "closer than the trees. Lower the merge distance or detect with a larger ",
      "window."), distance),
      call. = FALSE)
  }
  if (!is.null(chm)) {
    chm <- as.matrix(chm)
    h_tops <- .height_at(chm, tops)
    h_out <- .height_at(chm, out)
    lost <- is.na(h_out)
    # not above hmin: where detect_tops() could not have put a top
    if (!is.null(hmin)) lost <- lost | !(h_out > hmin)
    for (k in which(lost)) {
      idx <- which(g == k)              # rowsum(reorder = TRUE) keeps 1, 2, ...
      if (length(idx) > 1L && any(!is.na(h_tops[idx])))
        out[k, ] <- tops[idx[which.max(h_tops[idx])], ]
    }
  }
  out
}

# CHM value under each top; NA for a top outside the raster or on no data.
.height_at <- function(chm, tops) {
  rc <- .cells(tops)
  inside <- rc[, 1] >= 1 & rc[, 1] <= nrow(chm) & rc[, 2] >= 1 & rc[, 2] <= ncol(chm)
  inside[is.na(inside)] <- FALSE
  h <- rep(NA_real_, nrow(tops))
  if (any(inside)) h[inside] <- chm[cbind(rc[inside, 1], rc[inside, 2])]
  h
}

#' Drop tree tops standing on low CHM values
#'
#' Useful as a second, stricter pass: detect at 2 m so the growing can see
#' the understorey, then keep only tops above 10 m.
#'
#' @param chm Numeric matrix, the CHM the tops were detected on.
#' @param tops Two-column matrix of (row, col) positions.
#' @param hmin Numeric. Minimum height to keep.
#'
#' @return A two-column matrix. Tops outside the raster are dropped rather
#'   than wrapping around it, and so are tops standing on a cell with no
#'   data: a top [merge_tops()] placed between two crowns can land on one,
#'   and the trees it stood for are then left without a crown.
#' @export
screen_tops <- function(chm, tops, hmin) {
  chm <- as.matrix(chm)
  tops <- .as_tops(tops)
  if (!is.numeric(hmin) || length(hmin) != 1L || is.na(hmin))
    stop("hmin must be a single number", call. = FALSE)
  if (nrow(tops) == 0L) return(tops)
  rc <- .cells(tops)
  inside <- rc[, 1] >= 1 & rc[, 1] <= nrow(chm) & rc[, 2] >= 1 & rc[, 2] <= ncol(chm)
  inside[is.na(inside)] <- FALSE      # a top with an NA coordinate is nowhere
  keep <- rep(FALSE, nrow(tops))
  if (any(inside)) {
    v <- chm[cbind(rc[inside, 1], rc[inside, 2])]
    # NA >= hmin is NA, and an NA in a logical index returns a row of NA
    # rather than dropping the row: the top would stay, as (NA, NA), and
    # grow_crowns() would stop on it. A merged top is the centroid of its
    # group, so on a CHM with holes it can stand on one. numpy arrives at
    # the same answer by a different route, since NaN >= x is False.
    keep[inside] <- !is.na(v) & v >= hmin
  }
  tops[keep, , drop = FALSE]
}

#' Convert subpixel tops to integer pixel positions
#'
#' Each top goes to the cell it stands on, the nearest one: whole numbers are
#' the centres of the cells, and a half goes up.
#'
#' @param tops Two-column matrix of (row, col) positions.
#' @return An integer two-column matrix suitable for [grow_crowns()].
#' @export
as_pixels <- function(tops) {
  tops <- .as_tops(tops)
  matrix(as.integer(.cells(tops)), ncol = 2)
}
