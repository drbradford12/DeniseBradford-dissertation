# Chapter 2 figures, taken from the November 2026 proposal (second_attempt.qmd):
#   fig_beeswarm_vs_tiespread.png  -- proposal chunk fig-beeswarm-vs-tiespread
#   fig_categorical_fix_[a-d].png  -- proposal chunk fig-categorical-fix (mtcars)
# Usage: Rscript scripts/make_ch2_figures.R figures
args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args)) args[1] else "figures"
dir.create(out, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(patchwork)
  library(ggbeeswarm); library(ggpcp)
})

# ---- Beeswarm vs. tie_spread (18-observation toy data) -----------------------
raw_bs <- tibble(
  id    = 1:18,
  yearA = c(1990, 1991, 1993, 1995, 1996, 1997,
            1998, 1999, 2000, 2001, 2002, 2003,
            2003, 2004, 2004, 2005, 2006, 2007),
  yearB = c(1980, 1980, 1980, 1980, 1980, 1980,
            2000, 2000, 2000, 2000, 2000, 2000,
            2001, 2001, 2001, 2001, 2001, 2001),
  value = c(15,   40,   60,   25,   80,   50,
            35,   70,   20,   90,   45,   65,
            30,   55,   85,   10,   75,   95)
)
s01_bs <- function(x) (x - min(x)) / (max(x) - min(x))
dat_bs <- raw_bs %>% mutate(sA = s01_bs(yearA), sB = s01_bs(yearB), sV = s01_bs(value))

uB_bs     <- sort(unique(dat_bs$sB))
delta_bs  <- min(diff(uB_bs))
spread_bs <- delta_bs * 0.55

dat_broken_bs <- dat_bs %>%
  group_by(sB) %>%
  arrange(sA, .by_group = TRUE) %>%
  mutate(
    n    = n(),
    rank = row_number(),
    lo   = pmax(0, sB - spread_bs / 2),
    hi   = pmin(1, sB + spread_bs / 2),
    sBk  = if_else(n == 1L, sB, lo + (rank - 1L) / (n - 1L) * (hi - lo))
  ) %>%
  ungroup()

long_bs <- bind_rows(
  tibble(id = dat_broken_bs$id, x = 1, y = dat_broken_bs$sA,  sA = dat_broken_bs$sA),
  tibble(id = dat_broken_bs$id, x = 2, y = dat_broken_bs$sBk, sA = dat_broken_bs$sA),
  tibble(id = dat_broken_bs$id, x = 3, y = dat_broken_bs$sV,  sA = dat_broken_bs$sA)
)

p_bee <- ggplot(dat_bs, aes(x = factor(yearB), y = sA, colour = sA)) +
  geom_beeswarm(size = 3, cex = 3) +
  scale_colour_viridis_c(guide = "none") +
  labs(x = "Year B", y = "Year A (scaled)",
       subtitle = "ggbeeswarm: single axis, density-packed") +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        plot.subtitle = element_text(face = "bold"))

tie_boxes_bs <- dat_broken_bs %>%
  group_by(sB) %>%
  filter(n() > 1) %>%
  summarise(lo = first(lo), hi = first(hi), .groups = "drop") %>%
  transmute(xmin = 1.9, xmax = 2.1, ymin = lo, ymax = hi)

p_pcp_bs <- ggplot(long_bs, aes(x, y, group = id)) +
  geom_rect(data = tie_boxes_bs,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE,
            fill = "#AED6F1", colour = "#2471A3", alpha = 0.35, linewidth = 0.5) +
  geom_line(aes(colour = sA), linewidth = 0.6, show.legend = FALSE) +
  geom_point(aes(colour = sA), size = 1.4, show.legend = FALSE) +
  scale_colour_viridis_c() +
  scale_x_continuous(breaks = 1:3, labels = c("Year A", "Year B", "Value")) +
  labs(x = NULL, y = "scaled value",
       subtitle = "tie_spread: resolution-anchored + cross-axis sorted") +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        plot.subtitle = element_text(face = "bold"))

ggsave(file.path(out, "fig_beeswarm_vs_tiespread.png"),
       p_bee + p_pcp_bs + plot_layout(widths = c(1, 1.4)),
       width = 9, height = 4.5, dpi = 200)

# ---- mtcars: categorical tie-breaking, four panels ---------------------------
pcp_data <- mtcars %>%
  mutate(across(c(cyl, am, gear, carb), as.factor)) %>%
  pcp_select(cyl, am, gear, carb) %>%
  pcp_scale(method = "uniminmax")

panel <- function(d) {
  ggplot(d, aes_pcp()) +
    geom_pcp(aes(color = cyl), alpha = 0.6) +
    geom_pcp_boxes(boxwidth = 0.1, fill = "white", alpha = 0.5) +
    geom_pcp_labels() +
    theme_pcp() +
    guides(color = "none")
}

# (b) tie-breaking without sorting: shuffle positions within each level.
# The proposal used permute::shuffle() without a seed; a seed makes the
# panel reproducible.
set.seed(20261108)
pcp_shuffled <- mtcars %>%
  mutate(across(c(cyl, am, gear, carb), as.factor)) %>%
  pcp_select(cyl, am, gear, carb) %>%
  pcp_scale(method = "uniminmax") %>%
  pcp_arrange("from-right") %>%
  group_by(pcp_x, pcp_level) %>%
  mutate(order = sample(n()),
         pcp_y = pcp_y[order],
         pcp_yend = pcp_yend[order]) %>%
  ungroup()

panels <- list(
  a = panel(pcp_data),
  b = panel(pcp_shuffled),
  c = panel(pcp_data %>% pcp_arrange("from-left")),
  d = panel(pcp_data %>% pcp_arrange("from-right"))
)
for (k in names(panels)) {
  ggsave(file.path(out, sprintf("fig_categorical_fix_%s.png", k)), panels[[k]],
         width = 4.5, height = 3.2, dpi = 200)
}
cat("Chapter 2 figures written to", out, "\n")
