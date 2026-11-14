# =============================================================================
# Example: tie_spread numerical tie-breaking in the ggpcp pipeline
# (matches the July 2026 revision of the proposal)
#
# Run in your own R session (requires ggpcp, tidyverse, palmerpenguins,
# janitor, patchwork). Sources pcp_spread_ties.R from the same directory.
# =============================================================================
library(ggpcp)
library(tidyverse)
library(palmerpenguins)
library(janitor)
library(patchwork)

source(here::here("comprehensive_exam/numerical_tie_attempts/pcp_spread_ties.R"))
theme_set(theme_bw())

# Okabe-Ito qualitative palette, per the accessibility commitment in the
# proposal (Section 2, "Two design consequences follow for this work").
species_colors <- c("setosa" = "#E69F00", "versicolor" = "#CC79A7", "virginica" = "#009E73")
penguin_colors <- c("Adelie" = "#E69F00", "Chinstrap" = "#CC79A7", "Gentoo" = "#009E73")

# -----------------------------------------------------------------------------
# 1. IRIS -- the standard ggpcp pipeline with ONE extra step
# -----------------------------------------------------------------------------
iris_clean <- iris %>%
  na.omit() %>%
  clean_names()

iris_pcp <- iris_clean %>%
  pcp_select(sepal_length:species) %>%
  pcp_scale("uniminmax") %>%
  pcp_arrange("from-right") %>%                # categorical ties (ggpcp)
  pcp_spread_ties("from-right", sigma = 0.9)   # numerical ties (tie_spread)


# Annotation boxes (Section 4.5): light gray, anchored at [L(v), U(v)],
# labelled with the ORIGINAL tied values recovered from the raw data.
iris_boxes <- pcp_tie_boxes(iris_pcp, boxwidth = 0.2, original = iris_clean)

p_iris <- iris_pcp %>%
  ggplot(aes_pcp()) +
  geom_pcp_tieboxes(iris_boxes, labels = FALSE) +   # boxes BEHIND the lines
  geom_pcp_axes() +
  geom_pcp(aes(colour = species), alpha = 0.5) +
  geom_pcp_boxes(boxwidth = 0.2) +                  # categorical boxes
  geom_pcp_labels() +
  theme_pcp() +
  scale_color_manual("Species", values = species_colors) +
  labs(title = "Iris: resolution-anchored numerical tie-breaking",
       subtitle = expression(paste("Bands of height ", S(v) == delta[j],
                                   ", from-right hierarchical sorting"))) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

p_iris

# -----------------------------------------------------------------------------
# 2. PENGUINS -- one-call wrapper handling BOTH variable types (Section 4.6.4)
#    Integer-valued measurements create frequent ties; watch for the
#    saturation warning of Section 4.6.3.4 on heavily-tied axes.
# -----------------------------------------------------------------------------
penguins_clean <- penguins %>%
  na.omit() %>%
  clean_names()

peng_pcp <- penguins_clean %>%
  pcp_select(bill_length_mm:body_mass_g, island, species) %>%
  pcp_scale("uniminmax") %>%
  pcp_arrange_unified("from-right", sigma = 0.9)   # categorical + numerical

# min_n suppresses boxes for small bands to reduce clutter (Section 3.2)
peng_boxes <- pcp_tie_boxes(peng_pcp, boxwidth = 0.15,
                            original = penguins_clean, min_n = 3)

p_peng <- peng_pcp %>%
  ggplot(aes_pcp()) +
  geom_pcp_tieboxes(peng_boxes, labels = TRUE) +   # with original-value labels
  geom_pcp_axes() +
  geom_pcp(aes(colour = species), alpha = 0.5) +
  geom_pcp_boxes() +
  geom_pcp_labels() +
  theme_pcp() +
  scale_color_manual("species", values = penguin_colors) +
  labs(title = "Penguins: annotated numerical ties",
       subtitle = "Boxes suppressed below tie multiplicity 3 (min_n)") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

p_peng

# -----------------------------------------------------------------------------
# 3. Selective spreading (per-axis scoping via `vars`, contribution (2))
# -----------------------------------------------------------------------------
peng_selective <- penguins_clean %>%
  pcp_select(bill_length_mm:body_mass_g, species) %>%
  pcp_scale("uniminmax") %>%
  pcp_arrange("from-right") %>%
  pcp_spread_ties("from-right",
                  vars = c("flipper_length_mm", "body_mass_g"))

# -----------------------------------------------------------------------------
# 4. from-both with an explicit center axis (Equation 1). The sweep starts
#    at the center (ties there spread in input order, or seeded from one
#    nominated neighbor via center_reference) and proceeds outward, each
#    axis referencing the neighbor nearer the center.
# -----------------------------------------------------------------------------
iris_both <- iris_clean %>%
  pcp_select(sepal_length:species) %>%
  pcp_scale("uniminmax") %>%
  pcp_arrange("from-both") %>%
  pcp_spread_ties("from-both", center = "petal_length",
                  center_reference = "petal_width")

# -----------------------------------------------------------------------------
# 5. Before/after comparison panel (cf. Figure 13 of the proposal)
# -----------------------------------------------------------------------------
make_panel <- function(data, title) {
  data %>%
    ggplot(aes_pcp()) +
    geom_pcp_axes() +
    geom_pcp(aes(colour = species), alpha = 0.5) +
    geom_pcp_labels() +
    theme_pcp() +
    scale_color_manual(values = species_colors) +
    labs(title = title) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
}

before <- iris_clean %>%
  pcp_select(sepal_length:species) %>%
  pcp_scale("uniminmax") %>%
  pcp_arrange("from-right")

after <- before %>%
  pcp_spread_ties("from-right", sigma = 0.9)


print(make_panel(before, "Without numerical tie-breaking") |
      make_panel(after,  "With tie_spread (sigma = 0.9)"))

# -----------------------------------------------------------------------------
# 6. Exploring sigma ahead of the calibration study (Section 3.2). Values
#    approaching 1 shrink the inter-band gap below perceptibility and will
#    trigger the Section 2 warning via min_gap.
# -----------------------------------------------------------------------------
sigma_panels <- lapply(c(0.5, 0.7, 0.9), function(s) {
  before %>%
    pcp_spread_ties("from-right", sigma = s) %>%
    make_panel(title = bquote(sigma == .(s)))
})

wrap_plots(sigma_panels, nrow = 1)

