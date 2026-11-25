# Figure to compare corrected and noncorrected values

# Load packages ----
library(tidyverse)
library(patchwork)

# Read in plots ----
HFC <- readRDS("figures/panels/fig5d_HFC.rds")
LFC <- readRDS("figures/panels/fig5b_LFC.rds")

HFC_nc <- readRDS("results/experimental_validation/HFC_nc/HFC_nc_fig5d.rds")
LFC_nc <- readRDS("results/experimental_validation/LFC_nc/LFC_nc_fig5b.rds")

# Edits to plots
p1 <- LFC +
  theme(legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.box.spacing = unit(0, "pt"),
        legend.spacing.y = unit(0, "pt")) +
  guides(color = guide_legend(order = 1),
         shape = guide_legend(order = 1))

p2 <- LFC_nc +
  theme(legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.box.spacing = unit(0, "pt"),
        legend.spacing.y = unit(0, "pt")) +
  guides(color = guide_legend(order = 1),
         shape = guide_legend(order = 1))

p3 <- HFC +
  theme(legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.box.spacing = unit(0, "pt"),
        legend.spacing.y = unit(0, "pt")) +
  guides(color = guide_legend(order = 1),
         shape = guide_legend(order = 1))

p4 <- HFC_nc +
  theme(legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.box.spacing = unit(0, "pt"),
        legend.spacing.y = unit(0, "pt")) +
  guides(color = guide_legend(order = 1),
         shape = guide_legend(order = 1))

## Final plot ----
final_plot <- p1 + p2 + p3 + p4 + 
  plot_layout(
    ncol = 2,
    axes = "collect",
    axis_titles = "collect",
    guides = "collect"
  ) &
  theme(
    legend.position = "top"
  ) 

ggsave("figures/rfig_exp_comparison.pdf",
       final_plot, width = 7, height = 4.5, units = "in")
