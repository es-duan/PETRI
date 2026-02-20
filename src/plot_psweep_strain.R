# Plot parameter sweep data with relevant strain parameters on top. 

# Load packages ----
library(tidyverse)
library(argparse)
library(ggnewscale)
library(jsonlite)

# Set arguments parser inputs ----
parser <- ArgumentParser()
parser$add_argument("-p","--psweepsetting", type = "character", help = "Specify Treatment and parameter sweep setting")
parser$add_argument("-c","--colors", help = "JSON string of plot colors")

# Parse arguments
args <- parser$parse_args()

# Get treatment
ps <- args$psweepsetting
output_folder <- paste("results/parameter_sweeps", ps, sep = "/")

# Read in files ----
i1 <- readRDS(paste0(output_folder, "/", ps, "_inv_change_plot.rds"))

setting_list <- readRDS(paste0(output_folder, "/", ps, "_settings.rds"))
psi_ref <- as.numeric(setting_list$Ref_psi)

ph <- read_csv("input_data/strain_phenotypes.csv")

# Select relevant strains to plot ----
ps_s <- str_split_1(ps, "_")[2]
ps_sp <- str_sub(ps_s, 1,1)

ph_s <- ph %>%
  filter(str_detect(Strain, ps_sp)) %>%
  filter(Strain != ps_s) %>%
  mutate(log_conj = log10(Conjugation_rate))

## Plot strains on tile plot ----
i1_s <- i1 +
  new_scale_fill() +
  geom_point(data = ph_s,
             mapping = aes(log_conj, Growth_rate, fill = Strain),
             shape = 23, color = "white", size = 3)

# Save file
ggsave(paste0(output_folder, "/", ps, "_inv_change_strain_plot.pdf"),
       i1_s, height = 6.5, width = 8.775, units = "in")

