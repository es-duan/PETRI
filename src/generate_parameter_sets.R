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
total_points <- (Resolution)^2 #total points
print(paste("The total number of points for setting",ps,"is:",total_points))

## Determine values ----
log_gamma_range <- seq(log10(Ref_gamma)-Range_gamma, log10(Ref_gamma) + Range_gamma,
                       length.out = Resolution)
gamma_range <- 10^(log_gamma_range)
psi_range <- seq(Ref_psi - Range_psi, Ref_psi + Range_psi, 
                 length.out = Resolution)

## Pair values ----
# Pair gamma and psi fold change values
gamma_sweep <- rep(gamma_range,length(psi_range))
psi_sweep <- c(rep(psi_range[1],length(gamma_range)))

for(psi in 2:length(psi_range)){
  c = c(rep(psi_range[psi],length(gamma_range)))
  psi_sweep <- c(psi_sweep,c)
}

# Save final dataset
sweep_param <- data.frame(gamma_M = gamma_sweep,
                      psi_M = psi_sweep) %>%
  mutate(SetID = ps) %>%
  mutate(SweepID = row_number()) %>%
  select(SetID, SweepID, gamma_M, psi_M)

# Save files ----
saveRDS(setting_list, paste0(output_folder, "/", ps, "_settings.rds"))
write_csv(sweep_param, paste0(output_folder, "/", ps, "_params.csv"))
