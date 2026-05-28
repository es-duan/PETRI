# Plot parameter sweep data with relevant strain parameters on top. 

# Load packages ----
library(tidyverse)
library(argparse)
library(jsonlite)

# Set arguments parser inputs ----
parser <- ArgumentParser()
parser$add_argument("-p","--psweepsetting", type = "character", help = "Specify Treatment and parameter sweep setting")
parser$add_argument("-o","--points", help = "JSON string of point aesthetics")

# Parse arguments
args <- parser$parse_args()

# Load global variables ----
## Points ----
plot_points <- jsonlite::fromJSON(args$points)
psweep_point_size <- plot_points[["psweep_point_size"]]
sh_Anc <- plot_points[["sh_Anc"]]
sh_Mut <- plot_points[["sh_Mut"]]
sh_R1 <- plot_points[["sh_R1"]]
sh_copA <- plot_points[["sh_copA"]]
sh_finO <- plot_points[["sh_finO"]]

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
ps_sp <- str_split_1(ps_s, "-")[1]

ph_s <- ph %>%
  filter(str_detect(Strain, ps_sp)) %>%
  #filter(Strain != ps_s) %>%
  mutate(log_conj = log10(Conjugation_rate))

## Plot strains on tile plot ----
i1_s <- i1 +
  geom_point(data = ph_s,
             mapping = aes(log_conj, Growth_rate, shape = Strain),
             color = "black", size = psweep_point_size) +
  scale_shape_manual(values = c("S.pB10" = sh_Anc,
                                "S.pB10-A" = sh_Mut,
                                "S.pB10-B" = sh_finO,
                                "E.R1" = sh_R1,
                                "E.R1-copA" = sh_copA,
                                "E.R1-finO" = sh_finO))

# Save file
ggsave(paste0(output_folder, "/", ps, "_inv_change_strain_plot.pdf"),
       i1_s, height = 2.5, width = 3.75, units = "in")

