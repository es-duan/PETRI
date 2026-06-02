# Generate Fig Dimitriu psweep

# Load packages ----
library(tidyverse)
library(patchwork)
library(argparse)
library(jsonlite)

# Set arguments parser inputs ----
parser <- ArgumentParser()
parser$add_argument("-o","--points", help = "JSON string of point aesthetics")

# Parse arguments
args <- parser$parse_args()

# Load global variables ----
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
                                "E.R1-finO" = sh_finO),
                     labels = c("E.R1" = "wt",
                                "E.R1-copA" = expression(italic(copA)),
                                "E.R1-finO" = expression(italic(finO)^"-")),
                     name = "Genotype") +
  theme(legend.title = element_text(size = 8),
        legend.spacing.y = unit(0, "pt"),
        legend.spacing.x = unit(0, "pt")) +
  guides(shape = guide_legend(order = 1),
         fill = guide_colourbar(display = "gradient",
                                ticks = FALSE,
                                theme = theme(legend.text = element_text(size = 6)),
                                barwidth = unit(0.15, "in"), barheight = unit(0.6, "in")))

# Save plot ----
ggsave("figures/fig7_Dim_psweep.pdf",
       plot, width = 4, height = 3, units = "in")
