# Generate De Gelder mutant invasion figure

# Load packages ----
library(tidyverse)
library(patchwork)

# Read in figures ----
DG_dens <- readRDS("results/case_study_sims/DG_S.pB10/DG_S.pB10_density_plot.rds")
DG_freq <- readRDS("results/case_study_sims/DG_S.pB10/DG_S.pB10_frequency_plot.rds")

# Edits to plots ----
pA <- DG_dens +
  guides(linetype = "none") +
  # theme(legend.position = "none")
  theme(legend.position = "top",
        legend.title = element_blank(),
        legend.key.size = unit(0.15, "in"))

pB <- DG_freq +
  # theme(legend.position = "bottom",
  #       legend.title = element_blank())
  theme(legend.position = "none")


# Combine plots ----
final_plot <- pA + pB +
  plot_layout(nrow = 2,
              axes = "collect",
              axis_titles = "collect") +
  plot_annotation(tag_levels = "A")

ggsave("figures/fig4_DG_invasion.pdf",
       final_plot, width = 7.5, height = 5, units = "in")
