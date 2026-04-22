# Produce simrun case study plots

# Load packages ----
library(tidyverse)
library(argparse)
library(colorspace)

# Set arguments parser inputs ----
parser <- ArgumentParser()
parser$add_argument("-t","--treatment", type = "character", help = "Specify Treatment ID")

# Parse arguments
args <- parser$parse_args()

# Get treatment
treatment <- args$treatment

# Read in files ----
treatment_folder <- paste("results", "case_study_sims", treatment, sep = "/")

sim_dens <- read_csv(paste0(treatment_folder, "/", treatment, "_density_plot_df.csv"))
sim_freq <- read_csv(paste0(treatment_folder, "/", treatment, "_frequency_plot_df.csv"))
phases_rect <- read_csv(paste0(treatment_folder, "/", treatment, "_phases_plot_df.csv"))

# Set common aesthetics ----
# Colors
p_Anc <- "#8394F6"
p_Mut <- "#8A407A"
p_F <- "gray40"
p_C <- "gray80"
p_growth <- "#FFFBEF"
p_conj <- "#D0D9FF"
p_tselect <- "#FBE9FF"
p_imm <- "#d0ffeb"
p_select <- "#ffe9fc"

rifR_l <- "solid"
nalR_l <- "22"

# Retrieve ggplot theme
source("src/ggplot_theme.R")

# Plot densities over time ----
p1 <- ggplot() + 
  geom_rect(data = phases_rect,
            mapping = aes(xmin = T_start, xmax = T_end, ymin = 1, ymax = Inf, fill = Phase)) +
  scale_fill_manual(values = c("Growout" = "white",
                               "Growth" = p_batch,
                               "Conjugation" = p_conj,
                               "Transconjugant selection" = p_T)) +
  geom_line(data = sim_dens,
            mapping = aes(x = Time, y = Density,
                          color = Genotype_plot, linetype = Host),
            size = 2) +
  scale_color_manual(values = c("Ancestor" = p_Anc,
                                "Mutant" = p_Mut,
                                "Plasmid-free" = p_F,
                                "C" = p_C),
                     breaks = c("Ancestor", "Mutant", "Plasmid-free", "C"),
                     labels = c("Ancestor", "Mutant", "Plasmid-free", "C")) +
  scale_linetype_manual(values = c("rifR" = "solid",
                                   "nalR" = "dashed",
                                   "C" = "dotted")) +
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
                               "Growth" = p_batch,
                               "Conjugation" = p_conj,
                               "Transconjugant selection" = p_T)) +
  geom_line(data = sim_freq,
            mapping = aes(Time, Frequency, color = Genotype),
            size = 2) +
  scale_color_manual(values = c("Ancestor" = p_Anc,
                                "Mutant" = p_Mut)) +
  fig_aes

# Save plot
ggsave(paste0(treatment_folder, "/", treatment, "_frequency_plot.pdf"),
       p2, height = 5, width = 20, units = "in")
