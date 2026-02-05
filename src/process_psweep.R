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
t_s <- args$psweepsetting
sweep_folder <- paste("results/parameter_sweeps", t_s, sep = "/")

treatment <- str_extract(t_s, ".*?(?=_)")
setting <- str_extract(t_s, "(?<=_).*")

# Read in files ----
sweep_out <- read_csv(paste0(sweep_folder, "/" , t_s, "_out.csv"))
treatment_csv <- read.csv("input_data/Treatments_parameter_sweep.csv", header = F)

# Treatment file settings
tr_colnames <- treatment_csv[[1]]
treatment_csv <- as.data.frame(t(treatment_csv[,-1]))
colnames(treatment_csv) <- tr_colnames
row_number <- which(treatment_csv$Treatment_ID == treatment)

# Relevant Treatment file values
A1_0 = as.numeric(treatment_csv$A1_0[row_number])
M1_0 = as.numeric(treatment_csv$M1_0[row_number])

# Designate invasion status of mutations
sweep_plot <- sweep_out %>%
  mutate(Anc = A1 + A2,
         Mut = M1 + M2) %>%
  mutate(Anc_freq = Anc/(Anc + Mut),
         Mut_freq = Mut/(Mut + Anc)) %>%
  mutate(Mut_freq0 = M1_0/(A1_0 + M1_0)) %>%
  mutate(Mut_freq_change = ifelse(Mut_freq > Mut_freq0, "Increase", "Decrease")) %>%
  mutate(Invasion = case_when(Anc == 0 ~ "Displaced",
                              Mut== 0 ~ "Excluded",
                              # Designate situations where plasmid was outcompeted as NA
                              Anc == 0 & Mut == 0 ~ "NA",
                              Mut_freq > Mut_freq0 & Anc != 0 ~ "P_increase",
                              Mut_freq < Mut_freq0 & Anc != 0 ~ "P_decrease"))
# Save file
write_csv(sweep_plot, paste0(sweep_folder, "/", t_s, "_plot.csv"))
