# Process and plot strain phenotyping data

# Load packages ----
library(tidyverse)
library(broom)

# Read in data ----
gr_out <- read_csv("results/phenotyping/growth_rate/plating/growthrate0-4_av.csv")
LDM_av <- read_csv("results/phenotyping/LDM_conjugation/LDM_conjugation_av.csv")

# Rename data for merging ----
OD_growth <- gr_out %>%
  rename(Growth_rate_mean = gr_mean,
         Growth_rate_sd = gr_sd,
         Growth_rate_se = gr_se)

# Merge growth and conjugation data ----
phenotyping <- left_join(OD_growth, LDM_av, by = "Strain") %>%
  rename(Genotype = Strain)

# Export data ----
write_csv(phenotyping, "results/phenotyping/phenotyping_av.csv")

