# Generate full invasion simulation supplementary figure

# Load packages ----
library(tidyverse)
library(patchwork)

p_Anc <-  "#8A407A"
p_Mut <- "#8394F6"

# Read in figures ----
LFC_full <- readRDS("results/case_study_sims/LFC_S.pB10-A_full/LFC_S.pB10-A_full_frequency_plot.rds")
HFC_full <- readRDS("results/case_study_sims/HFC_S.pB10_full/HFC_S.pB10_full_frequency_plot.rds")

# Edits to plots ----
pA <- LFC_full +
  scale_color_manual(values = c(p_Anc, p_Mut),
                     labels = c("X", "Y")) +
  theme(legend.position = "top",
        legend.title = element_blank(),
        legend.key.size = unit(0.15, "in"))

pB <- HFC_full +
  # theme(legend.position = "bottom",
  #       legend.title = element_blank())
  theme(legend.position = "none")

# Combine plots ----
final_plot <- pA + pB +
  plot_layout(nrow = 2,
              axes = "collect",
              axis_titles = "collect") +
  plot_annotation(tag_levels = "A")

ggsave("figures/figS1_full_invasion.pdf",
       final_plot, width = 7, height = 5, units = "in")
