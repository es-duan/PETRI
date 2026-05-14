# Generate Fig Dimitriu psweep

# Load packages ----
library(tidyverse)
library(patchwork)
library(argparse)
library(jsonlite)

# Set arguments parser inputs ----
parser <- ArgumentParser()
parser$add_argument("-c","--colors", help = "JSON string of plot colors")
parser$add_argument("-o","--points", help = "JSON string of point aesthetics")

# Parse arguments
args <- parser$parse_args()

# Load global variables ----
## Colors ----
plot_colors <- jsonlite::fromJSON(args$colors)
p_imm <- plot_colors[["p_imm"]]

## Points ----
plot_points <- jsonlite::fromJSON(args$points)
psweep_point_size <- plot_points[["psweep_point_size"]]
sh_R1 <- plot_points[["sh_R1"]]
sh_copA <- plot_points[["sh_copA"]]
sh_finO <- plot_points[["sh_finO"]]

# Read in figures ----
Dim90_R1 <- readRDS("results/parameter_sweeps/pDim90_E.R1/pDim90_E.R1_inv_change_plot.rds")

# Strain data
ph <- read_csv("input_data/strain_phenotypes.csv")
ph_plot <- ph %>%
  filter(str_detect(Strain, "E")) %>%
  mutate(log_Conj = log10(Conjugation_rate))

# Modify plot ----
plot <- Dim90_R1 +
  scale_x_continuous(expand = c(0.025, 0.025)) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  geom_point(data = ph_plot,
             mapping = aes(log_Conj, Growth_rate, shape = Strain),
             size = psweep_point_size) +
  scale_shape_manual(values = c("E.R1" = sh_R1,
                                "E.R1-copA" = sh_copA,
                                "E.R1-finO" = sh_finO)) +
  guides(fill = guide_colourbar(theme = theme(legend.text = element_text(size = 10))))

# Save plot ----
ggsave("figures/fig7_Dim_psweep.pdf",
       plot, width = 7, height = 5, units = "in")
