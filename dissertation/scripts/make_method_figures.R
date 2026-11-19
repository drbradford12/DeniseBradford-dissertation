# Regenerate the Chapter 3 figures from the CURRENT pcp_spread_ties.R,
# using the code shown in numerical_ties_method.qmd.
# Usage: Rscript make_method_figures.R <path/to/pcp_spread_ties.R> <out_dir>
args <- commandArgs(trailingOnly = TRUE)
src <- args[1]; out <- args[2]
dir.create(out, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
  library(ggpcp); library(ggplot2); library(dplyr); library(patchwork)
})
source(src)

theme_set(theme_bw(base_size = 11))
species_colors <- c(setosa = "#E69F00", versicolor = "#CC79A7", virginica = "#009E73")

iris_clean <- iris
names(iris_clean) <- tolower(gsub("\\.", "_", names(iris_clean)))

iris_pcp <- iris_clean %>%
  pcp_select(sepal_length:species) %>%
  pcp_scale("uniminmax") %>%
  pcp_arrange("from-right") %>%
  pcp_spread_ties("from-right")

panel <- function(d, title) {
  ggplot(d, aes_pcp()) +
    geom_pcp_axes() +
    geom_pcp(aes(colour = species), alpha = 0.5, linewidth = 0.3) +
    geom_pcp_labels() +
    theme_pcp() +
    scale_color_manual("Species", values = species_colors) +
    labs(title = title, x = NULL, y = NULL) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1))
}

before <- iris_clean %>%
  pcp_select(sepal_length:species) %>%
  pcp_scale("uniminmax") %>%
  pcp_arrange("from-right")

after <- before %>% pcp_spread_ties("from-right", sigma = 0.9)

p1 <- panel(before, "(a) Without numerical tie-breaking") +
  panel(after, "(b) With tie-breaking (sigma = 0.9)")
ggsave(file.path(out, "fig_before_after.png"), p1, width = 10, height = 5, dpi = 200)

win_axes <- c("sepal_width", "petal_length", "petal_width")
win_y <- c((1.0 - 0.1) / 2.4, (1.5 - 0.1) / 2.4)   # petal widths 1.0 to 1.5 cm, scaled
win_pad <- 0.03
spread_win <- function(d, sigma) {
  suppressWarnings(pcp_spread_ties(d, "from-right", sigma = sigma))
}
win_panel <- function(d, title, extent = "inner", boxes_only = TRUE, labels = FALSE, min_n = 2) {
  dd <- d %>% filter(pcp_x %in% win_axes) %>% mutate(pcp_x = factor(pcp_x, levels = win_axes))
  bx <- pcp_tie_boxes(d, boxwidth = 0.3, min_n = min_n, extent = extent, original = iris_clean) %>%
    filter(pcp_x %in% win_axes) %>%
    mutate(axis_index = match(pcp_x, win_axes),
           xmin = axis_index - 0.15, xmax = axis_index + 0.15)
  g <- ggplot(dd, aes_pcp()) +
    { if (boxes_only) geom_pcp_tieboxes(bx) } +
    geom_pcp_axes() +
    geom_pcp(aes(colour = species), alpha = 0.5, linewidth = 0.3) +
    { if (labels) geom_pcp_tieboxes(bx, box = FALSE, labels = TRUE, label_size = 2.2,
                                    count_nudge = 0.17) } +
    scale_color_manual(values = species_colors, guide = "none") +
    coord_cartesian(xlim = c(0.75, 3.6), ylim = win_y + c(-win_pad, win_pad), expand = FALSE) +
    labs(title = title, x = NULL, y = NULL) +
    theme_bw(base_size = 10) +
    theme(panel.grid = element_blank(), axis.text.y = element_blank(),
          axis.ticks.y = element_blank(), plot.title = element_text(size = 10))
  g
}
full_panel <- function(d, extent = "inner") {
  dd <- d %>% filter(pcp_x %in% win_axes) %>% mutate(pcp_x = factor(pcp_x, levels = win_axes))
  bx <- pcp_tie_boxes(d, boxwidth = 0.3, extent = extent) %>%
    filter(pcp_x %in% win_axes) %>%
    mutate(axis_index = match(pcp_x, win_axes),
           xmin = axis_index - 0.15, xmax = axis_index + 0.15)
  ggplot(dd, aes_pcp()) +
    annotate("rect", xmin = 0.5, xmax = 3.5, ymin = win_y[1], ymax = win_y[2],
             fill = "#FFF3B0", alpha = 0.6) +
    geom_pcp_tieboxes(bx) +
    geom_pcp_axes() +
    geom_pcp(aes(colour = species), alpha = 0.5, linewidth = 0.3) +
    scale_color_manual(values = species_colors, guide = "none") +
    labs(x = NULL, y = NULL) +
    theme_bw(base_size = 10) +
    theme(panel.grid = element_blank(), axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
}
sigma_panels <- lapply(c(0.5, 0.7, 0.9), function(sg) {
  win_panel(spread_win(before, sg), sprintf("sigma = %.1f", sg))
})
p2 <- patchwork::wrap_plots(sigma_panels, nrow = 1)
ggsave(file.path(out, "fig_sigma.png"), p2, width = 10, height = 4, dpi = 200)

# box extent (proposal Fig. 15) and the two cues (proposal Fig. 16)
after9 <- spread_win(before, 0.9)
p_ext <- (win_panel(after9, "Box drawn to the allocation", extent = "band") +
            full_panel(after9, "band") + plot_layout(widths = c(1.6, 1))) /
         (win_panel(after9, "Box drawn to the dot-placement region", extent = "inner") +
            full_panel(after9, "inner") + plot_layout(widths = c(1.6, 1)))
ggsave(file.path(out, "fig_box_extent.png"), p_ext, width = 9, height = 7.5, dpi = 200)

p_cue <- (win_panel(after9, "Box with value and count labels", labels = TRUE, min_n = 5) +
            full_panel(after9) + plot_layout(widths = c(1.6, 1))) /
         (win_panel(after9, "Value and count labels alone", boxes_only = FALSE, labels = TRUE, min_n = 5) +
            full_panel(after9) + plot_layout(widths = c(1.6, 1)))
ggsave(file.path(out, "fig_two_cues.png"), p_cue, width = 9, height = 7.5, dpi = 200)
pw_delta <- min(diff(sort(unique((iris$Petal.Width - 0.1) / 2.4))))
cat("petal_width delta_j:", round(pw_delta, 3), " gaps at sigma .5/.7/.9:",
    round((1 - c(0.5, 0.7, 0.9)) * pw_delta, 3), "\n")

boxes <- pcp_tie_boxes(iris_pcp, boxwidth = 0.25, min_n = 5, original = iris_clean)
p3 <- ggplot(iris_pcp, aes_pcp()) +
  geom_pcp_tieboxes(boxes) +                                 # boxes BEHIND the lines
  geom_pcp_axes() +
  geom_pcp(aes(colour = species), alpha = 0.5, linewidth = 0.3) +
  geom_pcp_tieboxes(boxes, box = FALSE, labels = TRUE,
                    label_size = 2.2, count_nudge = 0.14) +  # labels in front
  geom_pcp_labels() +
  theme_pcp() +
  scale_color_manual("Species", values = species_colors) +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1))
ggsave(file.path(out, "fig_annotated.png"), p3, width = 8, height = 6, dpi = 200)

# Facts used in the chapter text and the repository guide
ties <- attr(iris_pcp, "pcp_ties")
cat("bands per axis:\n"); print(table(ties$pcp_x))
cat("band heights (ymax - ymin) by axis:\n")
print(tapply(ties$ymax - ties$ymin, ties$pcp_x, function(h) round(range(h), 4)))
cat("max multiplicity by axis:\n"); print(tapply(ties$n, ties$pcp_x, max))
cat("boxes drawn:", nrow(boxes), "\n")

# ---------------------------------------------------------------------------
# Palmer penguins application (layout of the November 2026 proposal, Figs.
# 16-17), drawn with pcp_spread_ties() rather than the proposal's inline
# re-implementation: 14 complete cases per species, from-left, sigma = 0.9.
# ---------------------------------------------------------------------------
suppressPackageStartupMessages(library(palmerpenguins))
axes <- c("bill_length_mm", "bill_depth_mm", "flipper_length_mm", "body_mass_g")
axis_labels <- c("Bill length\n(mm)", "Bill depth\n(mm)", "Flipper length\n(mm)", "Body mass\n(g)")
peng <- penguins %>%
  tidyr::drop_na(all_of(axes)) %>%
  group_by(species) %>% slice_head(n = 14) %>% ungroup()

peng_warnings <- character(0)
peng_pcp <- withCallingHandlers(
  peng %>%
    pcp_select(all_of(axes)) %>%
    pcp_scale("uniminmax") %>%
    ungroup() %>%                  # ungroup first: it would drop attr "pcp_ties"
    pcp_spread_ties("from-left", sigma = 0.9),
  warning = function(w) {
    peng_warnings <<- c(peng_warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  })

pb <- attr(peng_pcp, "pcp_ties")
species_cols <- c(Adelie = "#E69F00", Chinstrap = "#CC79A7", Gentoo = "#009E73")
half_w <- 0.05

peng_layers <- function() list(
  geom_pcp_axes(),
  geom_rect(data = pb, aes(xmin = axis_index - half_w, xmax = axis_index + half_w,
                           ymin = inner_min, ymax = inner_max),
            inherit.aes = FALSE, fill = "grey60", alpha = 0.6, colour = NA),
  geom_pcp(aes(colour = species), linewidth = 0.4, alpha = 0.9),
  scale_colour_manual(values = species_cols, name = NULL),
  scale_x_discrete(labels = setNames(axis_labels, axes)),
  labs(x = NULL, y = "scaled value"),
  theme_bw(base_size = 11),
  theme(panel.grid = element_blank(), legend.position = "top")
)

p_full <- ggplot(peng_pcp, aes_pcp()) + peng_layers()

# magnified window on the flipper-length axis around its largest band
zoom_axis <- 3
zb <- pb %>% filter(axis_index == zoom_axis) %>% slice_max(n, n = 1, with_ties = FALSE)
h  <- zb$inner_max - zb$inner_min
ylim <- c(zb$inner_min - 3 * h, zb$inner_max + 3 * h)
for (i in 1:5) {
  win <- pb %>% filter(axis_index == zoom_axis, inner_max >= ylim[1], inner_min <= ylim[2])
  new_ylim <- c(min(win$inner_min) - h, max(win$inner_max) + h)
  if (isTRUE(all.equal(new_ylim, ylim))) break
  ylim <- new_ylim
}
win_lab <- pcp_tie_boxes(peng_pcp, original = peng) %>%
  filter(axis_index == zoom_axis, inner_max >= ylim[1], inner_min <= ylim[2])

p_zoom <- ggplot(peng_pcp, aes_pcp()) + peng_layers() +
  geom_point(data = filter(peng_pcp, pcp_x == axes[zoom_axis]),
             aes(colour = species), size = 1.4) +
  coord_cartesian(xlim = c(zoom_axis - 0.3, zoom_axis + 0.3), ylim = ylim, expand = FALSE) +
  # only the zoomed axis gets a label; discrete breaks outside the window
  # would otherwise be drawn beside the panel
  scale_x_discrete(breaks = axes[zoom_axis], labels = axis_labels[zoom_axis]) +
  scale_y_continuous(breaks = (win_lab$inner_min + win_lab$inner_max) / 2,
                     labels = paste0(win_lab$label, " mm  (k = ", win_lab$n, ")")) +
  labs(y = NULL) +
  theme(legend.position = "none")

ggsave(file.path(out, "fig_penguins_spread.png"),
       p_full + p_zoom + plot_layout(widths = c(2.4, 1)),
       width = 9, height = 5, dpi = 200)

# value labels at the band centres, largest bands first, overlaps suppressed
lab <- pcp_tie_boxes(peng_pcp, original = peng) %>% arrange(desc(n))
p_lab <- p_full +
  geom_text(data = lab, aes(x = axis_index + half_w + 0.03, y = (ymin + ymax) / 2,
                            label = ifelse(n >= 3, paste0(label, " (n = ", n, ")"), label)),
            inherit.aes = FALSE, hjust = 0, size = 2.4, colour = "grey25",
            check_overlap = TRUE)
ggsave(file.path(out, "fig_penguins_labels.png"), p_lab, width = 9, height = 5, dpi = 200)

# facts quoted in the captions
scaled <- peng %>% select(all_of(axes)) %>%
  mutate(across(everything(), ~ (.x - min(.x)) / (max(.x) - min(.x))))
cat("penguins n =", nrow(peng), "\n")
cat("delta_j:", paste(names(scaled), round(sapply(scaled, function(y) min(diff(sort(unique(y))))), 4)), "\n")
cat("bands per axis:\n"); print(table(pb$pcp_x))
cat("max k per axis:\n"); print(tapply(pb$n, pb$pcp_x, max))
cat("zoom band:", zb$value, "k =", zb$n, "\n")
cat("warnings:\n"); cat(peng_warnings, sep = "\n")


# ---------------------------------------------------------------------------
# Flea beetles: the four treatments a user can apply (November 2026 proposal,
# Fig. 19), on the four many-valued measurements. Reconstructed from the
# figure and its caption; the numeric treatment uses pcp_spread_ties().
# ---------------------------------------------------------------------------
data(flea, package = "GGally")
fvars <- c("tars1", "tars2", "aede1", "aede3")
flea_cols <- c("Concinna " = "#56B4E9", "Heikert. " = "#E69F00", "Heptapot. " = "#009E73")
names(flea_cols) <- levels(flea$species)

flea_exact <- flea %>% pcp_select(all_of(fvars)) %>% pcp_scale("uniminmax") %>% ungroup()

# jitter at ggplot2's default: uniform offset within 40% of the axis resolution
set.seed(20261109)
flea_res <- sapply(fvars, function(a) {
  y <- flea[[a]]; min(diff(sort(unique(y)))) / diff(range(y))
})
flea_jit <- flea_exact %>%
  mutate(pcp_y = pcp_y + runif(n(), -0.4, 0.4) * flea_res[as.character(pcp_x)],
         pcp_yend = pcp_y)

flea_fac <- flea %>% mutate(across(all_of(fvars), as.factor)) %>%
  pcp_select(all_of(fvars)) %>% pcp_scale("uniminmax") %>% pcp_arrange("from-left") %>% ungroup()

flea_warnings <- character(0)
flea_ts <- withCallingHandlers(
  flea_exact %>% pcp_spread_ties("from-left", sigma = 0.9),
  warning = function(w) { flea_warnings <<- c(flea_warnings, conditionMessage(w)); invokeRestart("muffleWarning") })
fb <- attr(flea_ts, "pcp_ties")

# the beetles at the two largest gaps in the data
gap_info <- lapply(fvars, function(a) {
  u <- sort(unique(flea[[a]])); g <- diff(u); i <- which.max(g)
  side <- if (sum(flea[[a]] <= u[i]) < sum(flea[[a]] >= u[i + 1])) flea[[a]] <= u[i] else flea[[a]] >= u[i + 1]
  list(var = a, gap = max(g), ids = which(side))
})
gaps <- sapply(gap_info, `[[`, "gap"); names(gaps) <- fvars
top2 <- order(gaps, decreasing = TRUE)[1:2]
hi_ids <- unique(unlist(lapply(gap_info[top2], `[[`, "ids")))

flea_panel <- function(d, title, extra = NULL) {
  ggplot(d, aes_pcp()) + extra +
    geom_pcp_axes() +
    geom_pcp(aes(colour = species), linewidth = 0.3, alpha = 0.9) +
    geom_pcp(data = filter(d, pcp_id %in% hi_ids), colour = "black", linewidth = 1) +
    scale_colour_manual(values = flea_cols, guide = "none") +
    labs(title = title, x = NULL, y = NULL) +
    theme_bw(base_size = 10) +
    theme(panel.grid = element_blank(), axis.text.y = element_blank(),
          axis.ticks.y = element_blank(), plot.title = element_text(size = 10))
}
p_flea <- (flea_panel(flea_exact, "Exact values (no tie-breaking)") +
           flea_panel(flea_jit, "Jitter (ggplot2 default, 40% of resolution)")) /
          (flea_panel(flea_fac, "Factor coercion (pcp_arrange)",
                      geom_pcp_boxes(boxwidth = 0.1, fill = "white", alpha = 0.5)) +
           flea_panel(flea_ts, "pcp_spread_ties (variables stay numeric)",
                      geom_rect(data = fb, aes(xmin = axis_index - 0.05, xmax = axis_index + 0.05,
                                               ymin = inner_min, ymax = inner_max),
                                inherit.aes = FALSE, fill = "grey60", alpha = 0.7)))
ggsave(file.path(out, "fig_flea_treatments.png"), p_flea, width = 9, height = 8, dpi = 200)

cat("\nflea n =", nrow(flea), "\n")
cat("distinct values:", paste(fvars, sapply(flea[fvars], function(y) length(unique(y)))), "\n")
cat("beetles sharing a value:", paste(fvars, sapply(flea[fvars], function(y) sum(y %in% y[duplicated(y)]))), "\n")
cat("max multiplicity:", paste(fvars, sapply(flea[fvars], function(y) max(table(y)))), "\n")
cat("largest gaps:", paste(names(gaps), gaps), "| highlighted ids:", hi_ids, "\n")
cat("head distinct:", length(unique(flea$head)), " aede2 distinct:", length(unique(flea$aede2)), "\n")
cat("flea warnings:\n"); cat(flea_warnings, sep = "\n")


# ---------------------------------------------------------------------------
# Flea beetles, two regimes in one plot (November 2026 proposal, Fig. 23):
# aede2 (nine values) as a factor, the other five numeric, one unified call.
# ---------------------------------------------------------------------------
flea_two <- suppressWarnings(
  flea %>% mutate(aede2 = factor(aede2)) %>%
    pcp_select(species, 2:7, species) %>%
    pcp_scale("uniminmax") %>%
    pcp_arrange_unified("from-right", min_gap = 0))
flea_two_boxes <- pcp_tie_boxes(flea_two, boxwidth = 0.15, min_n = 4, original = flea)
p_two <- ggplot(flea_two, aes_pcp()) +
  geom_pcp_tieboxes(flea_two_boxes) +
  geom_pcp_axes() +
  geom_pcp(aes(colour = species), alpha = 0.7, linewidth = 0.35) +
  geom_pcp_tieboxes(flea_two_boxes, box = FALSE, labels = TRUE, label_size = 2.0,
                    count_nudge = 0.1) +
  geom_pcp_boxes(boxwidth = 0.15) +
  geom_pcp_labels() +
  scale_colour_manual(values = flea_cols, name = NULL) +
  theme_pcp() +
  theme(legend.position = "bottom")
ggsave(file.path(out, "fig_flea_regimes.png"), p_two, width = 10, height = 6.5, dpi = 200)
cat("two-regime boxes drawn:", nrow(flea_two_boxes), "\n")
