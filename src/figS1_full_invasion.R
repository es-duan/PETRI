# Generate full invasion simulation supplementary figure

# Load packages ----
library(tidyverse)
library(patchwork)
library(argparse)
library(jsonlite)

# Set arguments parser inputs ----
parser <- ArgumentParser()
parser$add_argument("-c","--colors", help = "JSON string of plot colors")

# Parse arguments
args <- parser$parse_args()

# Load global variables ----
## Colors ----
plot_colors <- fromJSON(args$colors)
p_Anc <- plot_colors[["p_Anc"]]
p_Mut <- plot_colors[["p_Mut"]]

## Retrieve ggplot theme ----
source("src/ggplot_theme.R")

# Read in figures ----
LFC_full <- readRDS("results/case_study_sims/LFC_S.pB10-A_full/LFC_S.pB10-A_full_frequency_plot.rds")
HFC_full <- readRDS("results/case_study_sims/HFC_S.pB10_full/HFC_S.pB10_full_frequency_plot.rds")

# Edits to plots ----
pA <- LFC_full +
  scale_color_manual(values = c(p_Anc, p_Mut),
                     labels = c("X", "Y")) +
  theme(legend.position = "top",
        legend.title = element_blank(),
        legend.key.size = unit(0.15, "in"),
        axis.title.x = element_blank())

pB <- HFC_full +
  theme(legend.position = "none")

# Combine plots ----
final_plot <- pA + pB +
  plot_layout(nrow = 2) +
  plot_annotation(tag_levels = "A")

ggsave("figures/figS1_full_invasion.pdf",
       final_plot, width = 7, height = 5, units = "in")
