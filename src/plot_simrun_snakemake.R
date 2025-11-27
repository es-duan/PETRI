# Produce simrun case study plots
# snakemake ver

# Load packages ----
library(tidyverse)
library(argparse)
library(jsonlite)

# Set arguments parser inputs ----
parser <- ArgumentParser()
parser$add_argument("-t","--treatment", type = "character", help = "Specify Treatment ID")
parser$add_argument("-c","--colors", help = "JSON string of plot colors")
parser$add_argument("-l","--lines", help = "JSON string of plot lines")

# Parse arguments
args <- parser$parse_args()

# Get treatment
treatment <- args$treatment

# Read in files ----
treatment_folder <- paste("results", "case_study_sims", treatment, sep = "/")

sim_dens <- read_csv(paste0(treatment_folder, "/", treatment, "_density_plot_df.csv"))
sim_freq <- read_csv(paste0(treatment_folder, "/", treatment, "_frequency_plot_df.csv"))
phases_rect <- read_csv(paste0(treatment_folder, "/", treatment, "_phases_plot_df.csv"))

# Load global variables ----
## Colors ----
plot_colors <- fromJSON(args$colors)
p_Anc <- plot_colors[["p_Anc"]]
p_Mut <- plot_colors[["p_Mut"]]
p_F <- plot_colors[["p_F"]]
p_growth <- plot_colors[["p_growth"]]
p_conj <- plot_colors[["p_conj"]]
p_tselect <- plot_colors[["p_tselect"]]
p_imm <- plot_colors[["p_imm"]]
p_select <- plot_colors[["p_select"]]

## Lines ----
plot_lines <- fromJSON(args$lines)
rifR_l <- plot_lines[["rifR_l"]]
nalR_l <- plot_lines[["nalR_l"]]

## Retrieve ggplot theme ----
source("src/ggplot_theme.R")

# Plot densities over time ----
p1 <- ggplot() + 
  geom_rect(data = phases_rect,
            mapping = aes(xmin = T_start, xmax = T_end, ymin = 1, ymax = Inf, fill = Phase)) +
  scale_fill_manual(values = c("Growout" = "white",
                               "Growth" = p_growth,
                               "Conjugation" = p_conj,
                               "Transconjugant selection" = p_tselect,
                               "Immigration" = p_imm,
                               "Plasmid selection" = p_select)) +
  geom_line(data = sim_dens,
            mapping = aes(x = Time, y = Density,
                          color = Genotype_plot, linetype = Host),
            size = 2) +
  scale_color_manual(values = c("Ancestor" = p_Anc,
                                "Mutant" = p_Mut,
                                "Plasmid-free" = p_F),
                     breaks = c("Ancestor", "Mutant", "Plasmid-free"),
                     labels = c("Ancestor", "Mutant", "Plasmid-free")) +
  scale_linetype_manual(values = c("rifR" = rifR_l,
                                   "nalR" = nalR_l)) +
  geom_hline(yintercept = 1.5, color = "white") +
  labs(x = "Time (hr)",
       y = "Density (CFU/mL)",
       color = "Genotype",
       fill = "Phase",
       linetype = "Host") +
  scale_y_continuous(trans = "log10",limits=c(1,2e9),
                     breaks=c(1e0,1e2,1e4,1e6,1e8),
                     labels=sapply(c(0,2,4,6,8),function(i){parse(text = sprintf("10^%d",i))}),
                     expand = c(0.01, 0.01)) +
  scale_x_continuous(expand = c(0.001, 0.001)) +
  guides(fill = guide_legend(order=1, reverse = TRUE),
         color = guide_legend(order=2),
         linetype = guide_legend(order=3)) +
  fig_aes

# Save plot
ggsave(paste0(treatment_folder, "/", treatment, "_density_plot.pdf"),
       p1, height = 5, width = 20, units = "in")

# Plot frequencies over time ----
p2 <- ggplot() + 
  geom_rect(data = phases_rect,
            mapping = aes(xmin = T_start, xmax = T_end, ymin = 0, ymax = Inf, fill = Phase)) +
  scale_fill_manual(values = c("Growout" = "white",
                               "Growth" = p_growth,
                               "Conjugation" = p_conj,
                               "Transconjugant selection" = p_tselect,
                               "Immigration" = p_imm,
                               "Plasmid selection" = p_select)) +
  geom_line(data = sim_freq,
            mapping = aes(Time, Frequency, color = Genotype),
            size = 2) +
  scale_color_manual(values = c("Ancestor" = p_Anc,
                                "Mutant" = p_Mut)) +
  scale_y_continuous(breaks=c(0, 0.000001,0.0001, 0.01, 1),
                     labels=c(0,sapply(c(-6,-4,-2),function(i){parse(text = sprintf("10^%d",i))}),1),
                     expand = c(0.01, 0.01),
                     trans = "log10") +
  scale_x_continuous(expand = c(0.001, 0.001)) +
  fig_aes

# Save plot
ggsave(paste0(treatment_folder, "/", treatment, "_frequency_plot.pdf"),
       p2, height = 5, width = 20, units = "in")
