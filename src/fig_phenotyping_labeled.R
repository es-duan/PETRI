# Annotate phenotyping plot

# Load packages ----
library(tidyverse)
library(ggpubr)
library(argparse)
library(jsonlite)

# Set arguments parser inputs ----
parser <- ArgumentParser()
parser$add_argument("-t","--treatment", type = "character", help = "Specify Treatment ID")
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

# Read in plot and phenotyping results ----
plot <- readRDS("results/phenotyping/pB10_phenotyping.rds")
phenotyping <- read_csv("results/phenotyping/phenotyping_av.csv")
growth_tt <- read_csv("results/phenotyping/growth_rate/growthratemax_tt.csv")
LDM_tt <- read_csv("results/phenotyping/LDM_conjugation/LDM_conjugation_tt.csv")

# Plot coordinates
anc_gr <- phenotyping$Growth_rate_mean[phenotyping$Genotype == "Anc"]
anc_conj <- phenotyping$Conjugation_rate_mean[phenotyping$Genotype == "Anc"]
anc_conj_se <- phenotyping$Conjugation_rate_se[phenotyping$Genotype == "Anc"]

mut_gr <- phenotyping$Growth_rate_mean[phenotyping$Genotype == "Mut"]
mut_gr_se <- phenotyping$Growth_rate_se[phenotyping$Genotype == "Mut"]
mut_conj <- phenotyping$Conjugation_rate_mean[phenotyping$Genotype == "Mut"]

p_growth <- growth_tt$p.value
p_conj <- LDM_tt$p.value

final <- plot +
  # Growth rate brackets
  annotate("segment", linewidth = 0.3, 
           x = anc_conj + anc_conj_se + 1e-12, xend = anc_conj + anc_conj_se + 1e-12,
           y = anc_gr, yend = mut_gr) +
  annotate("segment", linewidth = 0.3, 
           x = anc_conj + anc_conj_se + 5e-13, xend = anc_conj + anc_conj_se + 1e-12,
           y = anc_gr, yend = anc_gr) +
  annotate("segment", linewidth = 0.3, 
           x = anc_conj + anc_conj_se + 5e-13, xend = anc_conj + anc_conj_se + 1e-12,
           y = mut_gr, yend = mut_gr) +
  annotate("text",
           x = anc_conj + anc_conj_se + 2e-12, y = mean(c(anc_gr, mut_gr)),
           label = paste0("p = ", round(p_growth, 2)), angle = 270) +
  # Conjugation rate brackets
  annotate("segment", linewidth = 0.3, 
           x = anc_conj, xend = mut_conj,
           y = mut_gr + mut_gr_se + 0.02, yend = mut_gr + mut_gr_se + 0.02) +
  annotate("segment", linewidth = 0.3, 
           x = anc_conj, xend = anc_conj,
           y = mut_gr + mut_gr_se + 0.01, yend = mut_gr + mut_gr_se + 0.02) +
  annotate("segment", linewidth = 0.3, 
           x = mut_conj, xend = mut_conj,
           y = mut_gr + mut_gr_se + 0.01, yend = mut_gr + mut_gr_se + 0.02) +
  annotate("text",
           x = 10^mean(c(log10(anc_conj), log10(mut_conj))), y = mut_gr + mut_gr_se + 0.03,
           label = paste0("p = ", round(p_conj, 2)))

ggsave("figures/panels/fig3b_phenotyping.pdf",
       final, width = 7.5, height = 7, units = "in")
saveRDS(final, "figures/panels/fig3b_phenotyping.rds")



