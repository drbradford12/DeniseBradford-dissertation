# =============================================================================
# Tests for pcp_spread_ties() against the July 2026 paper revision.
# Base R only; builds ggpcp-style long data by hand.
# =============================================================================
source("/home/claude/work/pcp_spread_ties.R")

fails <- 0
check <- function(label, cond) {
  status <- if (isTRUE(cond)) "PASS" else { fails <<- fails + 1; "FAIL" }
  cat(sprintf("[%s] %s\n", status, label))
}

make_long <- function(W, classes = rep("numeric", ncol(W))) {
  axes <- colnames(W); n <- nrow(W)
  data.frame(
    pcp_id    = rep(seq_len(n), times = length(axes)),
    pcp_x     = factor(rep(axes, each = n), levels = axes),
    pcp_y     = as.vector(W),
    pcp_yend  = as.vector(W),
    pcp_class = rep(classes, each = n),
    stringsAsFactors = FALSE
  )
}
getY <- function(d, ax) {
  r <- as.character(d$pcp_x) == ax
  d$pcp_y[r][order(d$pcp_id[r])]
}

# -----------------------------------------------------------------------------
# Test 1: Figure 12 scenario. Variable B has tied groups at 0.8/0.5/0.2 with
# n = 5, 3, 2. from-right: initial axis is C (input order); B sorts by C.
# -----------------------------------------------------------------------------
set.seed(1)
n <- 12
B <- c(rep(0.8, 5), rep(0.5, 3), rep(0.2, 2), 0.95, 0.05)
A <- runif(n); C <- runif(n)          # A, C have no exact ties
long <- make_long(cbind(A = A, B = B, C = C))
res <- pcp_spread_ties(long, method = "from-right", sigma = 0.9)
Badj <- getY(res, "B")

delta <- 0.15; sigma <- 0.9           # distinct B values: .05 .2 .5 .8 .95

for (v in c(0.8, 0.5, 0.2)) {
  idx <- which(B == v); band <- Badj[idx]; k <- length(idx)
  # away from the boundary: S(v) = delta, inner region = sigma * delta,
  # centered at v (Eqs. 5-8)
  check(sprintf("v=%.1f band spans inner region [v - s*d/2, v + s*d/2]", v),
        isTRUE(all.equal(range(band),
                         c(v - sigma * delta / 2, v + sigma * delta / 2))))
  sp <- diff(sort(band))
  check(sprintf("v=%.1f uniform spacing Delta(v) = sigma*delta/(k-1)", v),
        isTRUE(all.equal(sp, rep(sigma * delta / (k - 1), k - 1))))
  # Step 4.b.iv: single reference neighbor = C (already processed, unchanged)
  check(sprintf("v=%.1f within-band rank matches reference axis C", v),
        all(order(band) == order(C[idx])))
}

check("untied values unchanged",
      isTRUE(all.equal(Badj[B %in% c(0.95, 0.05)], c(0.95, 0.05))))
check("axes without ties unchanged",
      isTRUE(all.equal(getY(res, "A"), A)) &&
        isTRUE(all.equal(getY(res, "C"), C)))
b08 <- range(Badj[B == 0.8]); b05 <- range(Badj[B == 0.5])
check("inter-band gap preserved (Gestalt separation)",
      b08[1] > b05[2] &&
        isTRUE(all.equal(b08[1] - b05[2], 0.3 - sigma * delta)))
check("pcp_yend mirrors pcp_y", isTRUE(all.equal(res$pcp_y, res$pcp_yend)))
tb <- attr(res, "pcp_ties")
check("pcp_ties metadata: 3 bands on B, allocation height = delta_j",
      nrow(tb) == 3 && all(tb$pcp_x == "B") &&
        isTRUE(all.equal(tb$ymax - tb$ymin, rep(delta, 3))))
check("pcp_ties records Delta(v) spacing",
      isTRUE(all.equal(tb$spacing, sigma * delta / (tb$n - 1))))

# -----------------------------------------------------------------------------
# Test 2: from-left references the LEFT neighbor A; axis A is the initial axis.
# -----------------------------------------------------------------------------
res_l <- pcp_spread_ties(long, method = "from-left", sigma = 0.9)
Badj_l <- getY(res_l, "B")
ok <- TRUE
for (v in c(0.8, 0.5, 0.2)) {
  idx <- which(B == v)
  ok <- ok && all(order(Badj_l[idx]) == order(A[idx]))
}
check("from-left: within-band rank matches reference axis A", ok)

# -----------------------------------------------------------------------------
# Test 3 (JULY SPEC): initial-axis convention + sequential keys.
# Ties on BOTH B and C. from-right: C is the INITIAL axis -> its band is
# spread in INPUT ORDER (Eq. 10), NOT by the left neighbor. B then sorts by
# C's ADJUSTED positions -- which are input-ordered, so B follows suit.
# -----------------------------------------------------------------------------
B2 <- c(rep(0.5, 4), 0.1, 0.9)
C2 <- c(rep(0.4, 4), 0.2, 0.8)
A2 <- c(0.9, 0.1, 0.7, 0.3, 0.5, 0.6)   # would give a DIFFERENT order than input
long2 <- make_long(cbind(A = A2, B = B2, C = C2))
res2 <- pcp_spread_ties(long2, method = "from-right", sigma = 0.9)
B2a <- getY(res2, "B"); C2a <- getY(res2, "C")
idx <- 1:4
check("initial axis (C): band spread in INPUT order, not by neighbor A",
      all(order(C2a[idx]) == 1:4) && !all(order(C2a[idx]) == order(A2[idx])))
check("sequential keys: B band ordered by ADJUSTED C positions",
      all(order(B2a[idx]) == order(C2a[idx])))

# Same data under from-left: A is initial (no ties there), C now references B.
res2l <- pcp_spread_ties(long2, method = "from-left", sigma = 0.9)
B2l <- getY(res2l, "B"); C2l <- getY(res2l, "C")
check("from-left: B references A; C references adjusted B",
      all(order(B2l[idx]) == order(A2[idx])) &&
        all(order(C2l[idx]) == order(B2l[idx])))

# -----------------------------------------------------------------------------
# Test 4: factor axes untouched; vars= restricts spreading.
# -----------------------------------------------------------------------------
Wm <- cbind(num1 = c(0.5, 0.5, 0.5, 1, 0), fac = c(.2, .4, .6, .8, 1),
            num2 = c(0.25, 0.25, 0.75, 0.75, 0.1))
lm_ <- make_long(Wm, classes = c("numeric", "factor", "numeric"))
res3 <- pcp_spread_ties(lm_, method = "from-right")
check("factor axis values untouched",
      isTRUE(all.equal(getY(res3, "fac"), unname(Wm[, "fac"]))))
res4 <- pcp_spread_ties(lm_, method = "from-right", vars = "num1")
check("vars=: only requested axis spread",
      !isTRUE(all.equal(getY(res4, "num1"), unname(Wm[, "num1"]))) &&
        isTRUE(all.equal(getY(res4, "num2"), unname(Wm[, "num2"]))))

# -----------------------------------------------------------------------------
# Test 5 (JULY SPEC): boundary clamping with PROPORTIONAL margins (Eq. 7).
# Tie at v = 0 with delta = 0.5: L = 0, U = 0.25, S = 0.25 < delta.
# Margins = (1-sigma)*S/2; inner region = sigma*S. No beta buffer.
# -----------------------------------------------------------------------------
Wb <- cbind(A = c(.1, .9, .5), B = c(0, 0, 0.5))   # distinct B: 0, 0.5 -> delta = .5
resb <- pcp_spread_ties(make_long(Wb), method = "from-right", sigma = 0.9)
Bb <- getY(resb, "B")
S_clamped <- 0.25
m <- (1 - 0.9) * S_clamped / 2
check("clamped band: realized allocation S(v) = U - L < delta_j",
      isTRUE(all.equal(range(Bb[1:2]), c(0 + m, S_clamped - m))))
tbb <- attr(resb, "pcp_ties")
check("clamped band metadata: [L, U] = [0, delta/2], spacing = sigma*S/(k-1)",
      isTRUE(all.equal(c(tbb$ymin, tbb$ymax), c(0, S_clamped))) &&
        isTRUE(all.equal(tbb$spacing, 0.9 * S_clamped / 1)))

# Single-unique-value axis: centered, warning
Ws <- cbind(A = c(.1, .9), B = c(0.5, 0.5))
warned <- FALSE
ress <- withCallingHandlers(
  pcp_spread_ties(make_long(Ws), method = "from-right"),
  warning = function(w) { warned <<- TRUE; invokeRestart("muffleWarning") })
check("single-unique-value axis: centered, warning issued",
      warned && isTRUE(all.equal(getY(ress, "B"), c(0.5, 0.5))))

# -----------------------------------------------------------------------------
# Test 6: tolerance grouping treats values within epsilon as tied.
# -----------------------------------------------------------------------------
Wt <- cbind(A = c(.1, .9, .5), B = c(0.3, 0.3 + 1e-14, 0.7))
rest <- pcp_spread_ties(make_long(Wt), method = "from-right", tol = 1e-10)
Bt <- getY(rest, "B")
check("near-identical values (within tol) grouped and spread",
      Bt[1] != Bt[2] && abs(mean(Bt[1:2]) - 0.3) < 1e-9)

# -----------------------------------------------------------------------------
# Test 7 (JULY SPEC): required warnings.
# (a) inter-band gap below min_gap: sigma = 0.999 makes (1-sigma)*delta tiny.
# (b) band saturation: many ties in a narrow band -> Delta(v) < min_spacing.
# -----------------------------------------------------------------------------
warns <- character(0)
tmp <- withCallingHandlers(
  pcp_spread_ties(long, method = "from-right", sigma = 0.999),
  warning = function(w) { warns <<- c(warns, conditionMessage(w))
  invokeRestart("muffleWarning") })
check("gap warning fires when (1-sigma)*delta_j < min_gap",
      any(grepl("Inter-band gap", warns)))

Wsat <- cbind(A = runif(60),
              B = c(rep(0.5, 50), 0.501, seq(0, 1, length.out = 9)))
warns2 <- character(0)
tmp <- withCallingHandlers(
  pcp_spread_ties(make_long(Wsat), method = "from-right", sigma = 0.9),
  warning = function(w) { warns2 <<- c(warns2, conditionMessage(w))
  invokeRestart("muffleWarning") })
check("saturation warning fires when Delta(v) < min_spacing",
      any(grepl("saturation", warns2)))

warns3 <- character(0)
tmp <- withCallingHandlers(
  pcp_spread_ties(long, method = "from-right", sigma = 0.9),
  warning = function(w) { warns3 <<- c(warns3, conditionMessage(w))
  invokeRestart("muffleWarning") })
check("no spurious warnings in the benign case", length(warns3) == 0)

# -----------------------------------------------------------------------------
# Test 8 (JULY SPEC): from-both with Equation 1 -- center axis in input order,
# outward sweeps reference the neighbor nearer the center; center= and
# center_reference= arguments.
# -----------------------------------------------------------------------------
# p = 3, center defaults to 2 (= B). B ties -> input order. A and C then
# reference B's adjusted positions.
B3 <- c(rep(0.5, 4), 0.1, 0.9)
A3 <- c(rep(0.3, 4), 0.6, 0.9)
C3 <- c(rep(0.4, 4), 0.2, 0.8)
long3 <- make_long(cbind(A = A3, B = B3, C = C3))
resfb <- pcp_spread_ties(long3, method = "from-both", sigma = 0.9)
A3a <- getY(resfb, "A"); B3a <- getY(resfb, "B"); C3a <- getY(resfb, "C")
check("from-both: center axis (B) band in INPUT order",
      all(order(B3a[idx]) == 1:4))
check("from-both: both sides reference the center's adjusted positions",
      all(order(A3a[idx]) == order(B3a[idx])) &&
        all(order(C3a[idx]) == order(B3a[idx])))

# center = "A" makes A initial (input order); B references A; C references B
resfb2 <- pcp_spread_ties(long3, method = "from-both", center = "A")
A3b <- getY(resfb2, "A"); B3b <- getY(resfb2, "B"); C3b <- getY(resfb2, "C")
check("from-both with center='A' reduces to from-left behavior (Eq. 1)",
      all(order(A3b[idx]) == 1:4) &&
        all(order(B3b[idx]) == order(A3b[idx])) &&
        all(order(C3b[idx]) == order(B3b[idx])))

# center_reference: seed center B's ordering from neighbor C's ORIGINAL
# positions (C processed after the center, so its pre-spread values apply)
resfb3 <- pcp_spread_ties(long3, method = "from-both",
                          center = "B", center_reference = "C")
B3c <- getY(resfb3, "B")
check("from-both: center_reference seeds center from nominated neighbor",
      all(order(B3c[idx]) == order(C3[idx])))

# -----------------------------------------------------------------------------
# Test 9: sigma scales band width linearly (away from boundaries).
# -----------------------------------------------------------------------------
res_s <- pcp_spread_ties(long, method = "from-right", sigma = 0.5)
w_small <- diff(range(getY(res_s, "B")[B == 0.8]))
w_big   <- diff(range(Badj[B == 0.8]))
check("sigma scales band width linearly",
      isTRUE(all.equal(w_small / w_big, 0.5 / 0.9)))

# -----------------------------------------------------------------------------
# Test 10: pcp_tie_boxes -- rect coords, original-value labels, min_n filter.
# -----------------------------------------------------------------------------
orig <- data.frame(A = round(A * 10),
                   B = c(rep(8, 5), rep(5, 3), rep(2, 2), 9.5, .5),
                   C = C)
boxes <- pcp_tie_boxes(res, boxwidth = 0.2, original = orig)
check("tie boxes: one row per band with rect coords + labels",
      nrow(boxes) == 3 &&
        all(c("xmin", "xmax", "ymin", "ymax", "label", "n") %in% names(boxes)) &&
        setequal(boxes$label, c("8", "5", "2")) &&
        isTRUE(all.equal((boxes$xmin + boxes$xmax) / 2, rep(2, 3))))
check("tie boxes: default extent is the inner region, height sigma * delta_j",
      isTRUE(all.equal(boxes$ymax - boxes$ymin, rep(sigma * delta, 3))) &&
        isTRUE(all.equal(min(diff(sort(c(boxes$ymin, boxes$ymax)))[c(2, 4)]),
                         0.3 - sigma * delta)))
boxes_b <- pcp_tie_boxes(res, boxwidth = 0.2, original = orig, extent = "band")
check("tie boxes: extent = 'band' gives the full allocation [L(v), U(v)]",
      isTRUE(all.equal(boxes_b$ymax - boxes_b$ymin, rep(delta, 3))) &&
        isTRUE(all.equal(boxes_b$ymin, tb$ymin)) && identical(boxes_b$label, boxes$label))
boxes_f <- pcp_tie_boxes(res, min_n = 3)
check("tie boxes: min_n suppresses low-multiplicity bands (Sec. 3.2)",
      nrow(boxes_f) == 2 && all(boxes_f$n >= 3))

# -----------------------------------------------------------------------------
# Test 10b: geom_pcp_tieboxes -- box, value label and count label are separate
# options; the count label defaults to on with the value labels, above
# count_min_n. Needs ggplot2; skipped otherwise.
# -----------------------------------------------------------------------------
if (requireNamespace("ggplot2", quietly = TRUE)) {
  n_layers <- function(...) length(geom_pcp_tieboxes(boxes, ...))
  check("tieboxes: default draws the boxes only", n_layers() == 1)
  check("tieboxes: labels = TRUE adds value labels and, by default, count labels",
        n_layers(labels = TRUE) == 3)
  check("tieboxes: counts = FALSE leaves the count label out",
        n_layers(labels = TRUE, counts = FALSE) == 2)
  check("tieboxes: counts = TRUE works without value labels",
        n_layers(counts = TRUE) == 2)
  check("tieboxes: box = FALSE returns the labels alone",
        n_layers(box = FALSE, labels = TRUE) == 2 && n_layers(box = FALSE) == 0)
  cnt <- geom_pcp_tieboxes(boxes, box = FALSE, counts = TRUE)[[1]]$data
  check("tieboxes: count label only for bands with n >= count_min_n (default 3)",
        nrow(cnt) == 2 && all(cnt$n >= 3))
  cnt2 <- geom_pcp_tieboxes(boxes, box = FALSE, counts = TRUE, count_min_n = 2)[[1]]$data
  check("tieboxes: count_min_n = 2 labels every band", nrow(cnt2) == 3)
  check("tieboxes: no count layer when no band reaches count_min_n",
        n_layers(counts = TRUE, count_min_n = 99) == 1)
} else {
  cat("[SKIP] geom_pcp_tieboxes tests (ggplot2 not installed)\n")
}

# -----------------------------------------------------------------------------
# Test 11: tie_spread alias exists and matches.
# -----------------------------------------------------------------------------
check("tie_spread() alias behaves identically",
      isTRUE(all.equal(
        tie_spread(long, method = "from-right", sigma = 0.9)$pcp_y,
        res$pcp_y)))

cat(sprintf("\n%s\n", if (fails == 0) "ALL TESTS PASSED" else
  paste(fails, "TEST(S) FAILED")))
if (fails > 0) quit(status = 1)
