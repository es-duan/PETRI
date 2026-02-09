# Process parameter sweep data and generate a preliminary plot

# Load packages ----
library(tidyverse)
library(argparse)

# Set arguments parser inputs ----
parser <- ArgumentParser()
parser$add_argument("-p","--psweepsetting", type = "character", help = "Specify Treatment and parameter sweep setting")

# Parse arguments
args <- parser$parse_args()

# Get treatment
ps <- args$psweepsetting
output_folder <- paste("results", "parameter_sweeps", ps, sep = "/")

# Read in files ----
sweep_out <- read_csv(paste0(output_folder, "/" , ps, "_out.csv"))
setting_list <- readRDS(paste0(output_folder, "/", ps, "_settings.rds"))

# Relevant Treatment file values
A1_0 = as.numeric(setting_list$A1_0)
M1_0 = as.numeric(setting_list$M1_0)

# Designate invasion status of mutations
sweep_plot <- sweep_out %>%
  mutate(Anc = A1 + A2,
         Mut = M1 + M2) %>%
  mutate(Anc_freq = Anc/(Anc + Mut),
         Mut_freq = Mut/(Mut + Anc)) %>%
  mutate(Mut_freq0 = M1_0/(A1_0 + M1_0)) %>%
  mutate(Mut_freq_change = Mut_freq - Mut_freq0,
         Mut_freq_inv = ifelse(Mut_freq > Mut_freq0, "Increase", "Decrease"))

# Save file
write_csv(sweep_plot, paste0(output_folder, "/", ps, "_plot.csv"))
