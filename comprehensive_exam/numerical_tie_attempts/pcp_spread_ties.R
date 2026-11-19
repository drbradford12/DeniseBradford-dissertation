#' ============================================================================
#' tie_spread: Resolution-Anchored Numerical Tie-Breaking for ggpcp
#' ============================================================================
#'
#' Implements the tie_spread algorithm from:
#'   "A Grammar of Graphics Approach to Managing
#'   Numerical Ties in Parallel Coordinate Plots" (Sections 4.3, 4.6, 4.8).
#'
#' The algorithm (Section 4.8, Algorithm Summary):
#'
#'   Step 1  Scale all variables to [0, 1]                (done by pcp_scale)
#'   Step 2  Per axis, resolution delta_j = min gap between distinct values
#'   Step 3  Processing order by method; each axis references the single
#'           adjacent ALREADY-PROCESSED axis A_j (Equation 1):
#'             from-left : axes 1..p, A_j = j - 1; axis 1 is the initial axis
#'             from-right: axes p..1, A_j = j + 1; axis p is the initial axis
#'             from-both : outward from a center axis c, A_j = the neighbor
#'                         nearer the center; c is the initial axis
#'           The INITIAL axis of the sweep has no A_j: its ties are spread
#'           in input order (kappa_i = pi(i), Equation 10) -- no fallback to
#'           the opposite direction.
#'   Step 4  Per tied value v (steps 4.b.i-vi):
#'             i.   Clamped outer bounds L(v) = max(0, v - delta_j/2),
#'                                       U(v) = min(1, v + delta_j/2)
#'             ii.  Realized allocation  S(v) = U(v) - L(v)
#'             iii. Margins S+(v) = S-(v) = (1 - sigma) * S(v) / 2
#'                  (fractions of the REALIZED allocation, Equation 7, so the
#'                  clamped boundary case needs no special formula)
#'             iv.  Order T_j(v) by kappa_i = ytilde_{i, A_j}; input order on
#'                  the initial axis
#'             v.   Spacing Delta(v) = sigma * S(v) / (k - 1)   (Equation 9)
#'             vi.  Positions ytilde = L(v) + S-(v) + (m - 1) * Delta(v)
#'   Step 5  Return adjusted data
#'
#' Warnings required by the paper:
#'   * Section 2:       warn when the resolution is too small for the
#'                      inter-band gap (1 - sigma) * delta_j to be perceptible
#'   * Section 4.6.3.4: warn when high tie multiplicity saturates a band,
#'                      i.e. Delta(v) falls below a legible spacing
#'
#' Pipeline usage:
#'
#'   data %>%
#'     pcp_select(...) %>%
#'     pcp_scale("uniminmax") %>%
#'     pcp_arrange("from-right") %>%      # categorical ties (existing ggpcp)
#'     pcp_spread_ties("from-right")     # numerical ties   (this file)
#'
#' or via pcp_arrange_unified(), which performs both. `tie_spread` is
#' provided as an alias for pcp_spread_ties, matching the paper's name.
#'
#' Only base R is used internally, so the module adds no dependencies beyond
#' what ggpcp already imports. Complexity is O(p * n log n) as stated in the
#' paper: one sweep, sorting-dominated per axis, no revisiting.
#' ============================================================================


# ---------------------------------------------------------------------------
# Internal: group values into ties under floating-point tolerance epsilon
# (Table 1 / Step 4.a). Returns an integer group id per element; NA -> NA.
# ---------------------------------------------------------------------------
tie_groups <- function(y, tol = 1e-10) {
  g <- rep(NA_integer_, length(y))
  ok <- which(!is.na(y))
  if (length(ok) == 0L) return(g)
  o <- ok[order(y[ok])]
  ys <- y[o]
  new_grp <- c(TRUE, diff(ys) > tol)
  g[o] <- cumsum(new_grp)
  g
}

# ---------------------------------------------------------------------------
# Internal: axis resolution delta_j (Step 2, Equation 3) -- the smallest
# difference between representative values of distinct tie groups.
# ---------------------------------------------------------------------------
axis_resolution <- function(rep_values) {
  if (length(rep_values) < 2L) return(NA_real_)
  min(diff(sort(rep_values)))
}

empty_ties <- function() {
  data.frame(value = numeric(0), n = integer(0),
             ymin = numeric(0), ymax = numeric(0),
             inner_min = numeric(0), inner_max = numeric(0),
             spacing = numeric(0), pcp_x = character(0),
             axis_index = integer(0), stringsAsFactors = FALSE)
}

empty_boxes <- function() {
  data.frame(xmin = numeric(0), xmax = numeric(0),
             ymin = numeric(0), ymax = numeric(0),
             label = character(0), n = integer(0), pcp_x = character(0),
             stringsAsFactors = FALSE)
}


# ---------------------------------------------------------------------------
# Internal: spread the tied observations on ONE axis (Steps 4.a-4.b.vi).
#
#   y      numeric vector of scaled values on this axis (one per obs)
#   key    numeric vector of realized positions on the single reference
#          neighbor A_j, or NULL on the initial axis of the sweep
#          (Equation 10: kappa_i = input order). Ties in the key, and NAs,
#          fall back to input order for a stable, well-defined arrangement.
#   sigma  separation parameter in (0, 1]        (Equation 7)
#   tol    tie tolerance epsilon
#
# Returns list(y     = adjusted values,
#              ties  = per-band metadata for annotation boxes/diagnostics,
#              delta = axis resolution delta_j)
# ---------------------------------------------------------------------------
spread_one_axis <- function(y, key = NULL, sigma = 0.9, tol = 1e-10) {
  out <- y
  grp <- tie_groups(y, tol = tol)
  ok  <- !is.na(grp)

  rep_vals <- tapply(y[ok], grp[ok], function(v) v[1])
  delta_j  <- axis_resolution(as.numeric(rep_vals))

  sizes    <- table(grp[ok])
  tied_ids <- as.integer(names(sizes)[sizes > 1L])

  ties <- empty_ties()[, 1:7]

  if (length(tied_ids) == 0L)
    return(list(y = out, ties = ties, delta = delta_j))

  # Degenerate axis: a single unique value, so delta_j is undefined and no
  # band can be anchored. Ties stay centered at the value.
  if (is.na(delta_j)) {
    warning("Axis has a single unique value; ties left centered (no spread).",
            call. = FALSE)
    return(list(y = out, ties = ties, delta = delta_j))
  }

  for (gid in tied_ids) {
    idx <- which(!is.na(grp) & grp == gid)
    v   <- rep_vals[[as.character(gid)]]
    k   <- length(idx)

    # clamped band, realized height, and margins as proportions of the band;
    # at the boundary everything contracts proportionally
    L <- max(0, v - delta_j / 2)
    U <- min(1, v + delta_j / 2)
    S <- U - L
    margin <- (1 - sigma) * S / 2
    L_inner <- L + margin
    U_inner <- U - margin

    # order within the band by the reference axis; input order otherwise and
    # as the stable tiebreak
    if (!is.null(key)) {
      ord <- order(key[idx], idx, na.last = TRUE)
    } else {
      ord <- order(idx)
    }

    # equispaced positions across the inner region, endpoints inclusive
    out[idx[ord]] <- seq(L_inner, U_inner, length.out = k)

    ties <- rbind(ties, data.frame(
      value     = v,
      n         = k,
      ymin      = L,
      ymax      = U,
      inner_min = L_inner,
      inner_max = U_inner,
      spacing   = sigma * S / (k - 1)
    ))
  }

  list(y = out, ties = ties, delta = delta_j)
}

# ---------------------------------------------------------------------------
# pcp_spread_ties(): pipeline step that breaks NUMERICAL ties.
#
# Operates on the long-format data produced by pcp_select() |> pcp_scale()
# (and, typically, pcp_arrange() for the categorical axes). Expects the
# standard ggpcp columns: pcp_id, pcp_x, pcp_y, pcp_yend, pcp_class.
#
#   method    "from-right", "from-left", or "from-both" (Section 4.8 step 3).
#   sigma     separation parameter in (0, 1]; fraction of each realized
#             allocation used as the dot-placement band, with (1 - sigma)
#             reserved as the inter-band gap (default 0.9, provisional per
#             RQ2/Section 3.2 pending the calibration study).
#   tol       tie tolerance epsilon for floating-point comparison.
#   vars      optional character vector of axis names to spread; default NULL
#             spreads every axis with pcp_class "numeric"/"integer". This is
#             the per-axis selectivity of contribution (2) in Section 3.1.
#   center    for method = "from-both": the center axis at which the sweep
#             begins, given as an axis name or index. Default: ceiling(p/2).
#             (Equation 1; the choice of center is a user-facing decision.)
#   center_reference  for method = "from-both": optionally seed the center
#             axis's within-band ordering from one nominated neighbor axis
#             (name or index) instead of input order (Section 4.2).
#   min_gap   perceptibility threshold for the inter-band gap: warn when
#             (1 - sigma) * delta_j < min_gap on a spread axis (Section 2;
#             the JND-grounded default is future empirical work -- 0.005 of
#             the unit axis is a provisional stand-in).
#   min_spacing  legibility threshold for within-band spacing: warn when
#             Delta(v) < min_spacing, the high-multiplicity saturation
#             regime of Section 4.6.3.4 (provisional default 0.001).
#
# Band metadata is attached as attr(result, "pcp_ties") and consumed by
# pcp_tie_boxes() / geom_pcp_tieboxes() below.
# ---------------------------------------------------------------------------
pcp_spread_ties <- function(data,
                            method = c("from-right", "from-left", "from-both"),
                            sigma = 0.9, tol = 1e-10, vars = NULL,
                            center = NULL, center_reference = NULL,
                            min_gap = 0.005, min_spacing = 0.001) {
  method <- match.arg(method)
  stopifnot(is.numeric(sigma), length(sigma) == 1, sigma > 0, sigma <= 1)

  required <- c("pcp_id", "pcp_x", "pcp_y", "pcp_class")
  missing_cols <- setdiff(required, names(data))
  if (length(missing_cols) > 0L) {
    stop("`data` is missing column(s): ", paste(missing_cols, collapse = ", "),
         ". Run pcp_select() |> pcp_scale() (and pcp_arrange()) first.",
         call. = FALSE)
  }

  # --- axis bookkeeping ------------------------------------------------------
  if (is.factor(data$pcp_x)) {
    axes <- levels(data$pcp_x)
  } else {
    axes <- unique(as.character(data$pcp_x))
  }
  p <- length(axes)
  ax_chr <- as.character(data$pcp_x)

  classes <- vapply(axes, function(a) {
    as.character(data$pcp_class[ax_chr == a][1])
  }, character(1))

  numeric_axes <- axes[classes %in% c("numeric", "integer")]
  vars <- rlang::enquo(vars)
  if (!rlang::quo_is_null(vars)) {
    # tidyselect against the axes: a one-row frame whose columns are the axes
    axes_df <- as.data.frame(stats::setNames(replicate(p, NA, simplify = FALSE),
                                             axes), check.names = FALSE)
    selected <- names(tidyselect::eval_select(vars, data = axes_df))
    numeric_axes <- intersect(numeric_axes, selected)
  }
  if (length(numeric_axes) == 0L) {
    attr(data, "pcp_ties") <- empty_ties()
    return(data)
  }

  as_axis_index <- function(x, what) {
    if (is.null(x)) return(NULL)
    if (is.character(x)) {
      i <- match(x, axes)
      if (is.na(i)) stop("`", what, "` axis not found: ", x, call. = FALSE)
      i
    } else {
      i <- as.integer(x)
      if (i < 1L || i > p)
        stop("`", what, "` index out of range: ", x, call. = FALSE)
      i
    }
  }

  # --- wide matrices of current realized positions (rows = observations) ----
  # Y holds pcp_y (the left-facing position of an axis), Yend holds pcp_yend
  # (the right-facing position). They differ only for factor axes arranged
  # with pcp_arrange(method = "from-both").
  ids <- unique(data$pcp_id)
  Y <- matrix(NA_real_, nrow = length(ids), ncol = p,
              dimnames = list(NULL, axes))
  row_of <- match(data$pcp_id, ids)
  col_of <- match(ax_chr, axes)
  Y[cbind(row_of, col_of)] <- data$pcp_y
  Yend <- Y
  if ("pcp_yend" %in% names(data)) Yend[cbind(row_of, col_of)] <- data$pcp_yend

  # --- processing order and reference neighbour A_j --------------------------
  # A_j is the one adjacent axis already processed on the sweep. The initial
  # axis has A_j = NA and is spread in input order.
  if (method == "from-left") {
    order_idx <- seq_len(p)                                  # initial axis: 1
    A <- c(NA_integer_, seq_len(p - 1))                      # A_j = j - 1
  } else if (method == "from-right") {
    order_idx <- rev(seq_len(p))                             # initial axis: p
    A <- ifelse(order_idx < p, order_idx + 1L, NA_integer_)  # A_j = j + 1
  } else { # from-both: outward from center c
    c_idx <- as_axis_index(center, "center")
    if (is.null(c_idx)) c_idx <- as.integer(ceiling(p / 2))
    lower <- if (c_idx > 1) (c_idx - 1):1 else integer(0)
    upper <- if (c_idx < p) (c_idx + 1):p else integer(0)
    inter <- integer(0)
    for (i in seq_len(max(length(lower), length(upper)))) {
      if (i <= length(upper)) inter <- c(inter, upper[i])
      if (i <= length(lower)) inter <- c(inter, lower[i])
    }
    order_idx <- c(c_idx, inter)
    A <- vapply(order_idx, function(j) {
      if (j == c_idx) NA_integer_                 # center: no reference
      else if (j > c_idx) j - 1L                  # neighbour nearer center
      else j + 1L
    }, integer(1))
    # optional seeding of the center from a nominated neighbour
    cr <- as_axis_index(center_reference, "center_reference")
    if (!is.null(cr)) {
      if (abs(cr - c_idx) != 1L)
        stop("`center_reference` must be an axis adjacent to the center.",
             call. = FALSE)
      A[1] <- cr
    }
  }

  # --- single sweep, sequential processing -----------------------------------
  all_ties <- list()
  gap_warn_axes <- character(0)
  sat_warn <- character(0)

  for (s in seq_along(order_idx)) {
    j <- order_idx[s]
    axis_name <- axes[j]
    if (!(axis_name %in% numeric_axes)) next

    # Segments between axes j-1 and j run from pcp_yend on axis j-1 to pcp_y
    # on axis j: a reference to the left uses its right-facing position (yend),
    # a reference to the right uses its left-facing position (y).
    key <- if (is.na(A[s])) NULL
    else if (A[s] < j) Yend[, A[s]]
    else Y[, A[s]]

    res <- spread_one_axis(Y[, j], key = key, sigma = sigma, tol = tol)
    Y[, j] <- res$y
    Yend[, j] <- res$y

    if (nrow(res$ties) > 0L) {
      if (!is.na(res$delta) && (1 - sigma) * res$delta < min_gap) {
        gap_warn_axes <- c(gap_warn_axes, axis_name)
      }
      sat <- res$ties$spacing < min_spacing
      if (any(sat)) {
        sat_warn <- c(sat_warn, sprintf(
          "%s (value%s %s)", axis_name, if (sum(sat) > 1) "s" else "",
          paste(format(round(res$ties$value[sat], 3)), collapse = ", ")))
      }
      res$ties$pcp_x <- axis_name
      res$ties$axis_index <- j
      all_ties[[axis_name]] <- res$ties
    }
  }

  if (length(gap_warn_axes) > 0L) {
    warning("Inter-band gap (1 - sigma) * delta_j is below the perceptible ",
            "threshold (", min_gap, ") on axis/axes: ",
            paste(unique(gap_warn_axes), collapse = ", "),
            ". Adjacent tie bands may visually merge; consider a smaller ",
            "sigma or fewer displayed values.", call. = FALSE)
  }
  if (length(sat_warn) > 0L) {
    warning("Within-band spacing Delta(v) is below the legibility threshold (",
            min_spacing, ") -- band saturation at high tie multiplicity -- ",
            "for: ", paste(sat_warn, collapse = "; "),
            ". Individual observations may overplot inside these bands.",
            call. = FALSE)
  }

  # --- write adjusted positions back into the long data ----------------------
  spread_cols <- which(axes %in% numeric_axes)
  sel <- col_of %in% spread_cols
  data$pcp_y[sel] <- Y[cbind(row_of[sel], col_of[sel])]
  if ("pcp_yend" %in% names(data)) {
    data$pcp_yend[sel] <- data$pcp_y[sel]
  }

  ties_df <- if (length(all_ties) > 0L) {
    do.call(rbind, c(all_ties, list(make.row.names = FALSE)))
  } else empty_ties()
  attr(data, "pcp_ties") <- ties_df
  data
}

# Alias matching the algorithm's name in the paper.
tie_spread <- pcp_spread_ties

# ---------------------------------------------------------------------------
# pcp_arrange_unified(): drop-in replacement for pcp_arrange() that handles
# categorical ties via ggpcp::pcp_arrange() and then numerical ties via
# pcp_spread_ties(), realizing the unified framework of Section 4.6.4.
# ---------------------------------------------------------------------------
pcp_arrange_unified <- function(data, method = "from-right", space = 0.05,
                                .by_group = TRUE, sigma = 0.9, tol = 1e-10,
                                vars = NULL, center = NULL,
                                center_reference = NULL,
                                min_gap = 0.005, min_spacing = 0.001) {
  has_factor <- any(data$pcp_class == "factor")
  if (has_factor) {
    data <- pcp_arrange(data, method = method, space = space,
                        .by_group = .by_group)
  }
  spread_method <- if (method %in% c("from-both", "from-both-keep"))
    "from-both" else method
  pcp_spread_ties(data, method = spread_method, sigma = sigma, tol = tol,
                  vars = {{ vars }}, center = center,
                  center_reference = center_reference,
                  min_gap = min_gap, min_spacing = min_spacing)
}

# ---------------------------------------------------------------------------
# pcp_tie_boxes(): annotation-box pipeline (Section 4.5, Figure 17).
#
# Consumes the metadata attached by pcp_spread_ties() and returns a data
# frame ready for geom_rect()/geom_label():
#   xmin, xmax, ymin, ymax  -- box anchored at [L(v), U(v)] (Section 4.5)
#   label                   -- original tied value when `original` (a wide
#                              data frame aligned with pcp_id) is supplied,
#                              otherwise the scaled tied value
#   n                       -- tie multiplicity, for min_n filtering
#
#   min_n  suppress boxes for bands with fewer than min_n observations, the
#          clutter-reduction default flagged in Section 3.2 (default 2, i.e.
#          every band; a value with a single observation never gets a box).
# ---------------------------------------------------------------------------
pcp_tie_boxes <- function(data, boxwidth = 0.2, original = NULL,
                          digits = 2, min_n = 2) {
  ties <- attr(data, "pcp_ties")
  if (is.null(ties) || nrow(ties) == 0L) return(empty_boxes())

  ties <- ties[ties$n >= min_n, , drop = FALSE]
  if (nrow(ties) == 0L) return(empty_boxes())

  ties$xmin <- ties$axis_index - boxwidth / 2
  ties$xmax <- ties$axis_index + boxwidth / 2

  ax_chr <- as.character(data$pcp_x)
  fmt <- function(x) format(round(x, digits))

  # a representative member of each band: any observation inside the inner region
  member_of <- function(r) {
    which(ax_chr == ties$pcp_x[r] &
            data$pcp_y >= ties$inner_min[r] - 1e-12 &
            data$pcp_y <= ties$inner_max[r] + 1e-12)
  }

  ties$label <- vapply(seq_len(nrow(ties)), function(r) {
    member <- member_of(r)
    if (length(member) == 0L) return(fmt(ties$value[r]))
    if (!is.null(original)) {
      ax <- ties$pcp_x[r]
      if (ax %in% names(original)) {
        id <- data$pcp_id[member[1]]
        return(fmt(as.numeric(original[[ax]][id])))
      }
    }
    if ("pcp_level" %in% names(data)) {
      lev <- suppressWarnings(as.numeric(as.character(data$pcp_level[member[1]])))
      if (!is.na(lev)) return(fmt(lev))
    }
    fmt(ties$value[r])
  }, character(1))

  rownames(ties) <- NULL
  ties
}

# ---------------------------------------------------------------------------
# geom_pcp_tieboxes(): convenience ggplot2 layers drawing the light gray
# bounding boxes (and optional value labels) behind the polylines, mirroring
# ggpcp's categorical geom_pcp_boxes() convention (Section 4.5). Add BEFORE
# geom_pcp() so the boxes sit behind the lines. (Whether this becomes a
# ggproto extension of the categorical box geom is the Phase 1 decision.)
# ---------------------------------------------------------------------------
geom_pcp_tieboxes <- function(boxes, fill = "grey85", colour = NA,
                              alpha = 0.9, labels = FALSE,
                              label_fill = "grey40", label_colour = "white",
                              label_size = 2.5) {
  layers <- list(
    ggplot2::geom_rect(
      data = boxes,
      ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax,
                   ymin = .data$ymin, ymax = .data$ymax),
      fill = fill, colour = colour, alpha = alpha, inherit.aes = FALSE
    )
  )

  if (labels && nrow(boxes) > 0L) {
    layers <- c(layers, list(
      ggplot2::geom_label(
        data = boxes,
        ggplot2::aes(x = (.data$xmin + .data$xmax) / 2,
                     y = (.data$ymin + .data$ymax) / 2,
                     label = .data$label),
        fill = label_fill, colour = label_colour, size = label_size,
        label.padding = ggplot2::unit(0.1, "lines"),
        label.r = ggplot2::unit(0.1, "lines"),
        inherit.aes = FALSE
      )
    ))
  }

  layers
}
