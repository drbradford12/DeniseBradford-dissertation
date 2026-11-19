# =============================================================================
# Example: the tie_spread algorithm, step by step, in the ggpcp pipeline
# (follows the mathematical framework of the November 2026 proposal,
#  Section 4.6: Steps 1-6, Equations 1-10, Propositions 1-2, Table 2)
#
# `tie_spread` names the algorithm; `pcp_spread_ties()` is the function.
#
# Run in your own R session (requires ggpcp, tidyverse, palmerpenguins,
# janitor, patchwork, ggforce; GGally for the flea data in Section 10).
# Sources pcp_spread_ties.R from the same directory.
#
#   1. One extra pipeline step            Section 4.2
#   2. Diagnostics                        Section 3.2
#   3. The framework, step by step        Section 4.6.3, Steps 1-6
#   4. Properties of the construction     Section 4.6.3.7, Propositions 1-2
#   5. The separation parameter sigma     Equations 6-8
#   6. Processing modes, initial axis     Sections 4.2, 4.6.3.5
#   7. Per-axis spreading with `vars`     Section 3.1, contribution (2)
#   8. The visual cue, label placement    Section 4.5
#   9. The unified call                   Section 4.6.4, Table 2
#  10. Three regimes                      Section 4.1
# =============================================================================
library(ggpcp)
library(tidyverse)
library(palmerpenguins)
library(janitor)
library(patchwork)
library(ggforce)      # facet_zoom()

source(here::here("comprehensive_exam/numerical_tie_attempts/pcp_spread_ties.R"))
theme_set(theme_bw())

# Okabe-Ito qualitative palette, per the accessibility commitment in the
# proposal (Section 2, "Two design consequences follow for this work").
species_colors <- c("setosa" = "#E69F00", "versicolor" = "#CC79A7", "virginica" = "#009E73")
penguin_colors <- c("Adelie" = "#E69F00", "Chinstrap" = "#CC79A7", "Gentoo" = "#009E73")
flea_colors    <- c("#E69F00", "#56B4E9", "#009E73")

sigma       <- 0.9      # separation parameter, provisional default (RQ2)
min_gap     <- 0.005    # threshold of the inter-band gap diagnostic
min_spacing <- 0.001    # threshold of the saturation diagnostic

iris_clean <- iris %>%
  na.omit() %>%
  clean_names()

penguins_clean <- penguins %>%
  na.omit() %>%
  clean_names()

iris_axes <- c("sepal_length", "sepal_width", "petal_length", "petal_width")

make_panel <- function(data, title, colors = species_colors) {
  data %>%
    ggplot(aes_pcp()) +
    geom_pcp_axes() +
    geom_pcp(aes(colour = species), alpha = 0.5) +
    geom_pcp_labels() +
    theme_pcp() +
    scale_color_manual(values = colors) +
    labs(title = title) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
}

# -----------------------------------------------------------------------------
# Helpers for the visual cue (used from Section 8 on)
# -----------------------------------------------------------------------------
# cue_plot(): one plot built from the three elements of the cue, each of which
# geom_pcp_tieboxes() draws or leaves out:
#   box     the gray box, always behind the lines
#   labels  the value labels, "behind" the lines, in "front" of them, or "none"
#   counts  the count labels, "n = k", drawn with the value labels
# geom_pcp_tieboxes() draws its layers in the order they are added, so boxes
# go in one call before geom_pcp() and, for labels in front, the labels go in
# a second call after it. `min_n` thins the labels: only bands with at least
# that many observations are labelled.
cue_plot <- function(data, boxes, labels = c("front", "behind", "none"), min_n = 2,
                     colors = species_colors, boxwidth = 0.2, box = TRUE,
                     counts = TRUE) {
  labels <- match.arg(labels)
  label_layers <- geom_pcp_tieboxes(boxes[boxes$n >= min_n, , drop = FALSE],
                                    box = FALSE, labels = TRUE,
                                    counts = counts, count_min_n = min_n)
  p <- ggplot(data, aes_pcp())
  if (box) p <- p + geom_pcp_tieboxes(boxes)              # boxes only
  if (labels == "behind") p <- p + label_layers
  p <- p +
    geom_pcp_axes() +
    geom_pcp(aes(colour = species), alpha = 0.5)
  if (labels == "front") p <- p + label_layers
  if (any(data$pcp_class == "factor")) {                  # categorical boxes and labels
    p <- p + geom_pcp_boxes(boxwidth = boxwidth) + geom_pcp_labels()
  }
  p +
    theme_pcp() +
    scale_color_manual(values = colors) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
}

# keep_axes(): restrict a spread layout and its tie boxes to two or three
# ADJACENT axes. Positions are not recomputed, so every line segment is the
# one drawn in the full plot. ggforce::facet_zoom() does not zoom a discrete
# scale, and the axis variable of a ggpcp plot is discrete; plotting only the
# axes of interest is how the zoom is narrowed to them.
keep_axes <- function(data, boxes, axes) {
  all_axes <- levels(data$pcp_x)
  stopifnot(all(axes %in% all_axes), all(diff(match(axes, all_axes)) == 1))
  shift <- match(axes[1], all_axes) - 1
  list(
    data  = data %>%
      ungroup() %>%
      filter(pcp_x %in% axes) %>%
      mutate(pcp_x = droplevels(pcp_x)),
    boxes = boxes %>%
      filter(pcp_x %in% axes) %>%
      mutate(axis_index = axis_index - shift, xmin = xmin - shift, xmax = xmax - shift)
  )
}

# -----------------------------------------------------------------------------
# 1. ONE EXTRA PIPELINE STEP (Section 4.2)
#    pcp_scale() performs the scaling of Equation 1. pcp_arrange() breaks the
#    ties on categorical axes. pcp_spread_ties() is the fourth step and breaks
#    the ties on numeric axes.
# -----------------------------------------------------------------------------
before <- iris_clean %>%
  pcp_select(sepal_length:species) %>%
  pcp_scale("uniminmax") %>%                   # Equation 1
  pcp_arrange("from-right")                    # categorical ties (ggpcp)

# This call warns. The warning is the first diagnostic; Section 2 explains it.
iris_pcp <- before %>%
  pcp_spread_ties("from-right", sigma = sigma) # numerical ties (tie_spread)

print(make_panel(before,   "Exact values: ties collapse to points") |
        make_panel(iris_pcp, "tie_spread (sigma = 0.9)"))

# -----------------------------------------------------------------------------
# 2. DIAGNOSTICS (Section 3.2)
#    Two warnings are built. Both follow from quantities the framework gives
#    in closed form, so both can be predicted before the function is called.
# -----------------------------------------------------------------------------
# r_j is the smallest gap between adjacent distinct values, in data units.
# ggplot2::resolution(x, zero = FALSE) returns the same number for a double
# vector; for an integer vector it returns 1 whatever the gaps are, so the
# gap is computed directly here, as pcp_spread_ties() computes it.
resolution_table <- function(raw, axes, sigma) {
  tibble(axis = axes) %>%
    mutate(
      r_j     = map_dbl(axis, ~ min(diff(sort(unique(raw[[.x]]))))),
      range   = map_dbl(axis, ~ diff(range(raw[[.x]]))),
      delta_j = r_j / range,                             # Equation 2
      gap     = (1 - sigma) * delta_j,                   # Proposition 1 (iv)
      k_star  = floor(sigma * delta_j / min_spacing) + 2 # smallest saturating k
    )
}

# (a) Inter-band gap. The gap between adjacent bands is (1 - sigma) * delta_j.
#     The warning fires when it falls below `min_gap`. At sigma = 0.9 that is
#     every axis with delta_j < 0.05, which includes all four iris axes.
iris_res <- resolution_table(iris_clean, iris_axes, sigma)
print(iris_res)
stopifnot(all(iris_res$gap < min_gap))

# (b) Saturation (Section 4.6.3.4). The band does not grow with k, so
#     Delta(v) = sigma * delta_j / (k - 1) falls below `min_spacing` once
#     k reaches k_star. Body mass is on a 25 g grid over 3600 g: k_star = 8.
peng_axes <- c("bill_length_mm", "bill_depth_mm", "flipper_length_mm", "body_mass_g")
peng_res  <- resolution_table(penguins_clean, peng_axes, sigma) %>%
  mutate(
    max_k       = map_int(axis, ~ max(table(penguins_clean[[.x]]))),
    n_saturated = map2_int(axis, k_star, ~ sum(table(penguins_clean[[.x]]) >= .y))
  )
print(peng_res)

# Both warnings fire here, and the saturation warning names the same values.
peng_warned <- penguins_clean %>%
  pcp_select(bill_length_mm:body_mass_g, island, species) %>%
  pcp_scale("uniminmax") %>%
  pcp_arrange_unified("from-right", sigma = sigma)

# From here on the diagnostics are acknowledged: `min_gap = 0` silences the
# first and `min_spacing = 0` the second.

# -----------------------------------------------------------------------------
# 3. THE FRAMEWORK, STEP BY STEP (Section 4.6.3)
#    pcp_spread_ties() records every band in attr(, "pcp_ties"). Its columns
#    are the quantities of Steps 2-4.
# -----------------------------------------------------------------------------
iris_pcp <- before %>%
  pcp_spread_ties("from-right", sigma = sigma, min_gap = 0)

framework <- attr(iris_pcp, "pcp_ties") %>%
  as_tibble() %>%
  transmute(
    axis   = pcp_x,
    v      = value,                   # tied value on the scaled axis
    k      = n,                       # |T_j(v)|              Step 2, Eq. 3
    L      = ymin,                    # L(v)                  Step 3, Eq. 4
    U      = ymax,                    # U(v)
    S      = ymax - ymin,             # S(v) = U(v) - L(v)    Eq. 5
    margin = inner_min - ymin,        # S+(v) = S-(v)         Eq. 6
    inner  = inner_max - inner_min,   # sigma * S(v)          Eq. 7
    Delta  = spacing                  # Delta(v)              Step 4, Eq. 8
  )
print(framework)

# Step 1 (Eq. 2). One resolution per axis. Away from the axis ends every band
# has height S(v) = delta_j, whatever its multiplicity.
delta_used <- framework %>%
  group_by(axis) %>%
  summarise(delta_j = max(S), bands = n(), k_min = min(k), k_max = max(k))
print(delta_used)
stopifnot(all.equal(delta_used$delta_j[match(iris_axes, delta_used$axis)],
                    iris_res$delta_j))

# Step 3 (Eqs. 4-7). The margins are a fraction of the realized allocation,
# so a band clamped at an axis end needs no separate formula.
stopifnot(
  all.equal(framework$margin, (1 - sigma) * framework$S / 2),   # Eq. 6
  all.equal(framework$inner,  sigma * framework$S)              # Eq. 7
)
clamped <- framework %>%
  inner_join(select(delta_used, axis, delta_j), by = "axis") %>%
  filter(S < delta_j - 1e-12)
print(clamped)                        # bands at v = 0 or v = 1: S(v) = delta_j / 2

# Step 4 (Eq. 8). Spacing shrinks as k grows; the band does not.
stopifnot(all.equal(framework$Delta, sigma * framework$S / (framework$k - 1)))

# Steps 5 and 6 (Eqs. 9-10) on one band. Under from-right the reference
# neighbor A_j of petal_width is species; kappa_i is the position there.
positions <- function(d) {
  d %>%
    ungroup() %>%
    as_tibble() %>%
    transmute(pcp_id, axis = as.character(pcp_x), y = pcp_y)
}
wide <- function(d) pivot_wider(positions(d), names_from = axis, values_from = y)

band <- framework %>%
  filter(axis == "petal_width", k == 5, S > 0.03) %>%
  slice(1)
members <- wide(before) %>%
  filter(abs(petal_width - band$v) < 1e-10) %>%
  select(pcp_id) %>%
  inner_join(wide(iris_pcp), by = "pcp_id") %>%
  transmute(pcp_id, kappa = species, y_adj = petal_width) %>%   # Eq. 9
  arrange(kappa) %>%
  mutate(m = row_number(),
         eq_10 = band$L + band$margin + (m - 1) * band$Delta)   # Eq. 10
print(members)
stopifnot(all.equal(members$y_adj, members$eq_10, check.attributes = FALSE))

# -----------------------------------------------------------------------------
# 4. PROPERTIES OF THE CONSTRUCTION (Section 4.6.3.7)
#    The two propositions, checked on the iris layout.
# -----------------------------------------------------------------------------
moved <- positions(before) %>%
  inner_join(rename(positions(iris_pcp), y_adj = y), by = c("pcp_id", "axis")) %>%
  filter(axis %in% iris_axes) %>%
  inner_join(select(iris_res, axis, delta_j), by = "axis")

# Proposition 1 (ii): no observation moves more than delta_j / 2, and no more
# than sigma * delta_j / 2 away from the axis ends.
shift <- moved %>%
  mutate(at_end = y < delta_j / 2 | y > 1 - delta_j / 2) %>%
  group_by(axis) %>%
  summarise(delta_j       = first(delta_j),
            max_shift     = max(abs(y_adj - y)),
            max_shift_mid = max(abs(y_adj - y)[!at_end]))
print(shift)
stopifnot(all(shift$max_shift     <= shift$delta_j / 2 + 1e-12),
          all(shift$max_shift_mid <= sigma * shift$delta_j / 2 + 1e-12))

# Proposition 1 (iii): away from the axis ends, a band is centered on its value.
stopifnot(with(anti_join(framework, clamped, by = c("axis", "v")),
               all.equal((L + U) / 2, v)))

# Proposition 1 (i) and (iv): order is preserved, and adjacent values stay
# apart. Away from the axis ends the gap is at least (1 - sigma) * delta_j.
# Next to a clamped band the upper margin is halved with the allocation, and
# the gap is 3/4 of that.
gaps <- moved %>%
  group_by(axis, delta_j, y) %>%
  summarise(lo = min(y_adj), hi = max(y_adj), .groups = "drop_last") %>%
  arrange(y, .by_group = TRUE) %>%
  mutate(gap    = lo - lag(hi),
         at_end = lag(y) < delta_j / 2 | y > 1 - delta_j / 2) %>%
  filter(!is.na(gap)) %>%
  summarise(min_gap_mid = min(gap[!at_end]),
            min_gap_end = if (any(at_end)) min(gap[at_end]) else NA_real_,
            .groups = "drop") %>%
  mutate(bound = (1 - sigma) * delta_j)
print(gaps)
stopifnot(all(gaps$min_gap_mid > 0),                                 # (i)
          all(gaps$min_gap_mid >= gaps$bound - 1e-12),               # (iv)
          all(gaps$min_gap_end >= 0.75 * gaps$bound - 1e-12, na.rm = TRUE))

# Proposition 2: no two lines cross inside a band. Within each band, the order
# on the axis must agree with the order on the reference neighbor.
band_inversions <- function(spread, exact, axis, reference) {
  wide(exact) %>%
    transmute(pcp_id, v = .data[[axis]]) %>%
    inner_join(wide(spread), by = "pcp_id") %>%
    group_by(v) %>%
    filter(n() > 1) %>%
    arrange(.data[[axis]], .by_group = TRUE) %>%
    summarise(inversions = sum(diff(.data[[reference]]) < 0), .groups = "drop") %>%
    summarise(bands = n(), inversions = sum(inversions))
}

axes_all   <- levels(before$pcp_x)
from_right <- tibble(axis = iris_axes, reference = axes_all[match(iris_axes, axes_all) + 1])
crossings  <- from_right %>%
  mutate(map2_dfr(axis, reference, ~ band_inversions(iris_pcp, before, .x, .y)))
print(crossings)
stopifnot(all(crossings$inversions == 0))

# The guarantee is tied to the reference neighbor. Bands ordered from the left
# are free of crossings on their left side and not on their right.
iris_left <- before %>%
  pcp_spread_ties("from-left", sigma = sigma, min_gap = 0)
print(from_right %>%
        mutate(map2_dfr(axis, reference, ~ band_inversions(iris_left, before, .x, .y))))

# -----------------------------------------------------------------------------
# 5. THE SEPARATION PARAMETER sigma (Equations 6-8)
#    sigma * delta_j is used for the observations; (1 - sigma) * delta_j is
#    left as the gap between adjacent bands. The choice of default is RQ2.
# -----------------------------------------------------------------------------
sigma_panels <- lapply(c(0.5, 0.7, 0.9), function(s) {
  before %>%
    pcp_spread_ties("from-right", sigma = s, min_gap = 0, min_spacing = 0) %>%
    make_panel(title = bquote(sigma == .(s)))
})
print(wrap_plots(sigma_panels, nrow = 1))

# The trade-off in numbers, on petal_length (delta_j = 0.017): a wider gap
# buys separation between bands and costs spacing inside them.
print(tibble(sigma = c(0.5, 0.7, 0.9)) %>%
        mutate(gap    = (1 - sigma) * iris_res$delta_j[3],
               k_star = floor(sigma * iris_res$delta_j[3] / min_spacing) + 2))

# -----------------------------------------------------------------------------
# 6. PROCESSING MODES AND THE INITIAL AXIS (Sections 4.2 and 4.6.3.5)
#    The construction is the same in every mode; only the sweep and the
#    reference neighbor A_j change. The initial axis of a sweep has no A_j,
#    so its ties are spread in input order.
# -----------------------------------------------------------------------------
# from-left: axes 1..p, A_j is the left neighbor (built above as iris_left).

# from-both: outward from `center`, each axis referencing the neighbor nearer
# the center. `center_reference` seeds the center from one adjacent axis
# instead of input order.
iris_both <- iris_clean %>%
  pcp_select(sepal_length:species) %>%
  pcp_scale("uniminmax") %>%
  pcp_arrange("from-both") %>%
  pcp_spread_ties("from-both", center = "petal_length",
                  center_reference = "petal_width", min_gap = 0)

print(make_panel(iris_left, "from-left") | make_panel(iris_both, "from-both, center = petal_length"))

# -----------------------------------------------------------------------------
# 7. PER-AXIS SPREADING WITH `vars` (Section 3.1, contribution (2))
#    `vars` takes a tidyselect specification. Axes not selected, and factor
#    axes, are left at their values.
# -----------------------------------------------------------------------------
peng_selective <- penguins_clean %>%
  pcp_select(bill_length_mm:body_mass_g, species) %>%
  pcp_scale("uniminmax") %>%
  pcp_arrange("from-right") %>%
  pcp_spread_ties("from-right", vars = c(flipper_length_mm, body_mass_g),
                  min_gap = 0, min_spacing = 0)

unique(attr(peng_selective, "pcp_ties")$pcp_x)   # only the two selected axes

# -----------------------------------------------------------------------------
# 8. THE VISUAL CUE AND LABEL PLACEMENT (Section 4.5)
#    Three elements: the box says where the tie is and how far the spread
#    extends, the value label says what the tied value is, and the count label
#    says how many observations share it.
# -----------------------------------------------------------------------------
# A box covers the inner region of its band, the sigma * delta_j in which the
# observations sit, and is labelled with the ORIGINAL tied value. Every band
# is kept here (min_n = 2); labels are thinned with `min_n` in cue_plot()
# (Section 3.2).
iris_boxes <- pcp_tie_boxes(iris_pcp, boxwidth = 0.2, original = iris_clean)

# At full scale a band is a few percent of the axis at most. The zoom is set
# on three adjacent axes and on a range of values named on one of them: here
# petal widths 1.0 to 1.5 cm, where every value on the 0.1 cm grid is tied.
iris_axes3  <- c("sepal_width", "petal_length", "petal_width")
iris_zoom   <- keep_axes(iris_pcp, iris_boxes, iris_axes3)
iris_window <- iris_zoom$boxes %>%
  filter(pcp_x == "petal_width", between(as.numeric(label), 1.0, 1.5)) %>%
  with(range(ymin, ymax) + c(-0.01, 0.01))

# (a) The extent of the box. With extent = "band" the box is the full
#     allocation [L(v), U(v)]. Bands at adjacent values share an edge, so on
#     petal width the boxes merge into one gray column although the values
#     differ. With the default, extent = "inner", adjacent boxes are separated
#     by the gap (1 - sigma) * delta_j and distinct values read as distinct.
iris_band <- keep_axes(iris_pcp,
                       pcp_tie_boxes(iris_pcp, boxwidth = 0.2, original = iris_clean,
                                     extent = "band"),
                       iris_axes3)

p_band <- cue_plot(iris_band$data, iris_band$boxes, labels = "none") +
  labs(title = 'extent = "band": boxes at adjacent values merge') +
  facet_zoom(ylim = iris_window, zoom.size = 1.5)

p_inner <- cue_plot(iris_zoom$data, iris_zoom$boxes, labels = "none") +
  labs(title = 'extent = "inner" (default): one box per tied value') +
  facet_zoom(ylim = iris_window, zoom.size = 1.5)

print(p_band / p_inner)

# (b) Where the labels are drawn is a design decision, and part of RQ2.
#   Behind the lines: the lines stay traceable through the band, and a label
#     is hidden wherever many lines pass over it.
#   In front of the lines: the label is always legible, and it covers the
#     lines inside the band, which is where tracing happens.
# Both placements, on the same axes and in the same window:
p_behind <- cue_plot(iris_zoom$data, iris_zoom$boxes, labels = "behind", min_n = 5) +
  labs(title = "Labels behind the lines") +
  facet_zoom(ylim = iris_window, zoom.size = 1.5)

p_front <- cue_plot(iris_zoom$data, iris_zoom$boxes, labels = "front", min_n = 5) +
  labs(title = "Labels in front of the lines") +
  facet_zoom(ylim = iris_window, zoom.size = 1.5)

print(p_behind / p_front)

# (c) Whether to draw the box at all. The box is the provisional default, and
#     Section 4.5 sets out a perceptual case against it on a numeric axis:
#       - enclosure groups strongly, so tied values read as bins and the
#         metric axis as a set of levels (common region);
#       - a categorical box's height is a frequency, so a numeric box of
#         constant height invites an equal-count reading it cannot support;
#       - boxes add filled area and edges where the lines converge (clutter);
#       - their edges cross every line through a tie and may slow tracing.
#     With the value and count labels present, what the box alone still shows
#     is the extent of the spread. Both cues, in the same window, for RQ2:
p_box <- cue_plot(iris_zoom$data, iris_zoom$boxes, labels = "front", min_n = 5) +
  labs(title = "Box with value and count labels") +
  facet_zoom(ylim = iris_window, zoom.size = 1.5)

p_nobox <- cue_plot(iris_zoom$data, iris_zoom$boxes, labels = "front", min_n = 5,
                    box = FALSE) +
  labs(title = "Value and count labels alone") +
  facet_zoom(ylim = iris_window, zoom.size = 1.5)

print(p_box / p_nobox)

# (d) The count label is an option of geom_pcp_tieboxes(). By default it is
#     drawn whenever the value labels are (counts = labels), for boxes holding
#     at least `count_min_n` observations (3 unless changed):
#
#       geom_pcp_tieboxes(boxes, labels = TRUE)                   # values and counts
#       geom_pcp_tieboxes(boxes, labels = TRUE, counts = FALSE)   # values only
#       geom_pcp_tieboxes(boxes, counts = TRUE)                   # counts only
#       geom_pcp_tieboxes(boxes, labels = TRUE, count_min_n = 10) # counts for n >= 10
#
#     The band has the same height at every value, so without the count label
#     frequency can be read only from how densely a band is filled. With and
#     without it, in the same window:
p_counts <- cue_plot(iris_zoom$data, iris_zoom$boxes, labels = "front", min_n = 5) +
  labs(title = "counts = TRUE (default with value labels)") +
  facet_zoom(ylim = iris_window, zoom.size = 1.5)

p_nocounts <- cue_plot(iris_zoom$data, iris_zoom$boxes, labels = "front", min_n = 5,
                       counts = FALSE) +
  labs(title = "counts = FALSE") +
  facet_zoom(ylim = iris_window, zoom.size = 1.5)

print(p_counts / p_nocounts)

# The same call on the full plot, straight from the package function:
print(iris_pcp %>%
        ggplot(aes_pcp()) +
        geom_pcp_tieboxes(iris_boxes) +
        geom_pcp_axes() +
        geom_pcp(aes(colour = species), alpha = 0.5) +
        geom_pcp_tieboxes(filter(iris_boxes, n >= 10), box = FALSE, labels = TRUE) +
        geom_pcp_boxes(boxwidth = 0.2) +
        geom_pcp_labels() +
        theme_pcp() +
        scale_color_manual(values = species_colors) +
        labs(title = "Iris: boxes behind the lines, value and count labels in front",
             subtitle = "Labels shown for tied values shared by at least 10 flowers") +
        theme(legend.position = "none",
              axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)))

# -----------------------------------------------------------------------------
# 9. THE UNIFIED CALL (Section 4.6.4, Table 2)
#    One formula governs both axis types: spacing = (allocation - margins) /
#    (n - 1). On a categorical axis the allocation is set by frequency; on a
#    numeric axis by the resolution. pcp_arrange_unified() applies both.
# -----------------------------------------------------------------------------
peng_pcp <- penguins_clean %>%
  pcp_select(bill_length_mm:body_mass_g, island, species) %>%
  pcp_scale("uniminmax") %>%
  pcp_arrange_unified("from-right", sigma = sigma,   # categorical + numerical
                      min_gap = 0, min_spacing = 0)

peng_boxes <- pcp_tie_boxes(peng_pcp, boxwidth = 0.15, original = penguins_clean)

# Body mass is the saturation case of Section 2: the band at 3800 g holds 12
# birds in under 1% of the axis. The count label is the primary reading there,
# so the labels are drawn in front. The zoom takes the three adjacent
# measurements and body masses from 3700 to 3900 g.
peng_zoom   <- keep_axes(peng_pcp, peng_boxes,
                         c("bill_depth_mm", "flipper_length_mm", "body_mass_g"))
peng_window <- peng_zoom$boxes %>%
  filter(pcp_x == "body_mass_g", between(as.numeric(label), 3700, 3900)) %>%
  with(range(ymin, ymax) + c(-0.005, 0.005))

p_peng <- cue_plot(peng_zoom$data, peng_zoom$boxes, labels = "front", min_n = 8,
                   colors = penguin_colors, boxwidth = 0.15) +
  labs(title = "Penguins: three adjacent measurements, zoomed on 3700 to 3900 g",
       subtitle = "Labels suppressed below tie multiplicity 8 (min_n)") +
  facet_zoom(ylim = peng_window, zoom.size = 1.5)

print(p_peng)

# -----------------------------------------------------------------------------
# 10. THREE REGIMES (Section 4.1)
#     Which treatment a numeric variable needs depends on how many distinct
#     values it takes.
# -----------------------------------------------------------------------------
# Regime 3, near-continuous: nothing to spread. The axis is returned unchanged
# and no bands are recorded.
set.seed(2026)
continuous <- tibble(a = rnorm(50), b = rnorm(50)) %>%
  pcp_select(a, b) %>%
  pcp_scale("uniminmax")
continuous_spread <- pcp_spread_ties(continuous)
stopifnot(identical(continuous$pcp_y, continuous_spread$pcp_y),
          nrow(attr(continuous_spread, "pcp_ties")) == 0)

if (requireNamespace("GGally", quietly = TRUE)) {
  data(flea, package = "GGally")

  # aede2 takes 9 values and head 13; the other four take 33 to 46.
  print(summarise(flea, across(where(is.numeric), n_distinct)))

  # Regime 1, few integer levels: coerce to a factor. aede2 is treated that way.
  # Regime 2, many values with scattered ties: tie_spread. head sits at the
  # boundary and is kept numeric here.
  flea_regimes <- flea %>%
    mutate(aede2 = factor(aede2)) %>%
    pcp_select(species, tars1:aede3, species) %>%
    pcp_scale("uniminmax") %>%
    pcp_arrange_unified("from-right", sigma = sigma, min_gap = 0)

  # What coercion does in regime 2: every axis becomes a ladder of thin,
  # equally spaced boxes and the metric scale is gone.
  flea_coerced <- flea %>%
    mutate(across(tars1:aede3, factor)) %>%
    pcp_select(species, tars1:aede3, species) %>%
    pcp_scale("uniminmax") %>%
    pcp_arrange("from-right")

  flea_boxes <- pcp_tie_boxes(flea_regimes, boxwidth = 0.15, original = flea)

  flea_panel <- function(data, title, boxes = NULL) {
    p <- ggplot(data, aes_pcp())
    if (!is.null(boxes)) p <- p + geom_pcp_tieboxes(boxes, labels = FALSE)
    p +
      geom_pcp_axes() +
      geom_pcp(aes(colour = species), alpha = 0.7) +
      geom_pcp_boxes(boxwidth = 0.15) +
      theme_pcp() +
      scale_color_manual(values = flea_colors) +
      labs(title = title) +
      theme(legend.position = "none")
  }

  print(flea_panel(flea_coerced, "Every axis coerced to a factor") /
          flea_panel(flea_regimes, "aede2 as a factor; the other axes spread", flea_boxes))

  # The cue on three adjacent axes, zoomed on head widths 49 to 52. tars1 and
  # tars2 are in the many-valued regime; head, with 13 values, is the boundary.
  flea_zoom   <- keep_axes(flea_regimes, flea_boxes, c("tars1", "tars2", "head"))
  flea_window <- flea_zoom$boxes %>%
    filter(pcp_x == "head", between(as.numeric(label), 49, 52)) %>%
    with(range(ymin, ymax) + c(-0.01, 0.01))

  names(flea_colors) <- levels(flea$species)
  p_flea <- cue_plot(flea_zoom$data, flea_zoom$boxes, labels = "front", min_n = 3,
                     colors = flea_colors, boxwidth = 0.15) +
    labs(title = "Flea beetles: tars1, tars2 and head, zoomed on head widths 49 to 52") +
    facet_zoom(ylim = flea_window, zoom.size = 1.5)
  print(p_flea)
}
