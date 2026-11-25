# Script for sitching together Dim psweep figures

# Load packages ----
library(tidyverse)
library(patchwork)

# Read in figures ----
im0 <- readRDS("figures/panels/Dim_psweep0.rds")
im10 <- readRDS("figures/panels/Dim_psweep10.rds")
im68 <- readRDS("figures/panels/Dim_psweep68.rds")
im90 <- readRDS("figures/panels/Dim_psweep90.rds")

# Edits to plots ----
# text variables
tx <- -14.5
ty <- 1.1
tsize <- 3.5

pA <- im0 +
  annotate("text", x = tx, y = ty,
           label = "0% immigration",
           color = "white", hjust = 0, size = tsize) +
  guides(fill = guide_colourbar(position = "top",
                                display = "gradient",
                                ticks = FALSE,
                                theme = theme(legend.text = element_text(size = 6))),
         shape = "none")

pB <- im10 +
  annotate("text", x = tx, y = ty,
           label = "10% immigration",
           color = "white", hjust = 0, size = tsize) +
  guides(fill = guide_colourbar(position = "top",
                                display = "gradient",
                                ticks = FALSE,
                                theme = theme(legend.text = element_text(size = 6))),
         shape = "none")

pC <- im68 +
  annotate("text", x = tx, y = ty,
           label = "68% immigration",
           color = "white", hjust = 0, size = tsize) +
  guides(fill = guide_colourbar(position = "top",
                                display = "gradient",
                                ticks = FALSE,
                                theme = theme(legend.text = element_text(size = 6))),
         shape = "none")

pD <- im90 +
  annotate("text", x = tx, y = ty,
           label = "90% immigration",
           color = "white", hjust = 0, size = tsize) +
  guides(fill = guide_colourbar(position = "top",
                                display = "gradient",
                                ticks = FALSE,
                                theme = theme(legend.text = element_text(size = 6))),
         shape = guide_legend(position = "right"))

# Combine plots ----
final_plot <- pA + pB + pC + pD +
  plot_layout(nrow = 1,
              axes = "collect",
              axis_titles = "collect") +
  plot_annotation(tag_levels = "A") 

ggsave("figures/fig_Dim_psweeps.pdf",
       final_plot, width = 12, height = 4, units = "in")