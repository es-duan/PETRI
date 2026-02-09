# Script to generate the growth and conjugation parameters for a parameter sweep

# Load packages ----
library(tidyverse)
library(argparse)

# Set arguments parser inputs ----
parser <- ArgumentParser()
parser$add_argument("-p","--psweepsetting", type = "character", help = "Specify parameter sweep setting")

# Parse arguments
args <- parser$parse_args()
ps <- args$psweepsetting

# Load setting file ----
setting_csv <- read_csv(paste("input_data", "parameter_sweeps", 
                              paste0(ps, "_psweep_settings.csv"),
                              sep = "/"))

## Transpose file into a list for easier reading ----
s_colnames <- setting_csv[[1]]
setting_csv <- as.data.frame(t(setting_csv[,-1]))
colnames(setting_csv) <- s_colnames
setting_list <- as.list(setting_csv)

## Create folder to store results ----
output_folder <- paste("results", "parameter_sweeps", ps, sep = "/")

if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
} else{
  print(paste("Setting",ps,"folder exists. Rewriting previous run."))
}

# Set values for parameters ----
Ref_gamma  = as.numeric(setting_list$Ref_gamma)
Ref_psi  = as.numeric(setting_list$Ref_psi)

Range_gamma  = as.numeric(setting_list$Range_gamma)
Range_psi  = as.numeric(setting_list$Range_psi)

Resolution  = as.numeric(setting_list$Resolution)

# Calculate the set of mutant parameters to run ----

## Print the total number of values ----
total_points <- (2*Resolution)^2 #total points
print(paste("The total number of points for setting",ps,"is:",total_points))

## Determine fold change values ----
gamma_res <- Range_gamma/Resolution
gamma_range <- seq(-Range_gamma, Range_gamma, gamma_res)

psi_res <- Range_psi/Resolution
psi_range <- seq(-Range_psi, Range_psi, psi_res)

## Calculate fold change values ----
# Pair gamma and psi fold change values
gamma_sweep <- rep(gamma_range,length(psi_range))
psi_sweep <- c(rep(psi_range[1],length(gamma_range)))

for(psi in 2:length(psi_range)){
  c = c(rep(psi_range[psi],length(gamma_range)))
  psi_sweep <- c(psi_sweep,c)
}

fc_values <- data.frame(gamma_fold = gamma_sweep,
                      psi_change = psi_sweep)

# Calculate actual values
# Conjugation change is currently being calculated as a order of magnitude fold change
# Growth change is currently being calculated as an absolute change
sweep_param <- fc_values %>%
  mutate(gamma_M = Ref_gamma*(10^gamma_fold),
         psi_M = psi_change + Ref_psi) %>%
  mutate(SetID = ps) %>%
  mutate(SweepID = row_number()) %>%
  select(SetID, SweepID, gamma_fold, psi_change, gamma_M, psi_M)

# Save files ----
saveRDS(setting_list, paste0(output_folder, "/", ps, "_settings.rds"))
write_csv(sweep_param, paste0(output_folder, "/", ps, "_params.csv"))
