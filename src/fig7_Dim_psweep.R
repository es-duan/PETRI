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
parser$add_argument("-l","--lines", help = "JSON string of plot lines")

# Parse arguments
args <- parser$parse_args()

# Load global variables ----
## Colors ----
plot_colors <- jsonlite::fromJSON(args$colors)
p_Exc <- plot_colors[["p_Exc"]]
p_Dis <- plot_colors[["p_Dis"]]
p_lowI <- plot_colors[["p_lowI"]]
p_highI <- plot_colors[["p_highI"]]
p_lowD <- plot_colors[["p_lowD"]]
p_highD <- plot_colors[["p_highD"]]
p_axes <- plot_colors[["p_axes"]]
p_mid <- plot_colors[["p_mid"]]
p_invline <- plot_colors[["p_invline"]]

## Lines ----
plot_lines <- fromJSON(args$lines)
inv_width <- plot_lines[["inv_width"]]

## Points ----
plot_points <- jsonlite::fromJSON(args$points)
psweep_point_size <- plot_points[["psweep_point_size"]]
sh_R1 <- plot_points[["sh_R1"]]
sh_copA <- plot_points[["sh_copA"]]
sh_finO <- plot_points[["sh_finO"]]

## Retrieve ggplot theme ----
source("src/ggplot_theme.R")

# Read in data ----
sweep_plot <- read_csv("results/parameter_sweeps/pDim90_E.R1/pDim90_E.R1_plot.csv")

# Strain data
ph <- read_csv("input_data/strain_phenotypes.csv")
ph_plot <- ph %>%
  filter(str_detect(Strain, "E")) %>%
  mutate(log_Conj = log10(Conjugation_rate))

gamma_ref <- ph_plot$log_Conj[ph_plot$Strain == "E.R1"]
psi_ref <- ph_plot$Growth_rate[ph_plot$Strain == "E.R1"]

max_gamma <- gamma_ref + 4
min_gamma <- gamma_ref - 4

max_psi <- psi_ref + 0.4
min_psi <- psi_ref - 0.4

## Filter plots ----
sweep_plot2 <- sweep_plot %>% 
  filter(psi_M > min_psi & psi_M < max_psi) %>%
  filter(log_gamma_M > min_gamma & log_gamma_M < max_gamma)

# Re-calculate relative scaling
log_max_change <- max(sweep_plot2$log_Mut_freq_change)
log_min_change <- min(sweep_plot2$log_Mut_freq_change)

sweep_plot_f <- sweep_plot2 %>%
  mutate(log_Mut_freq_change2 = ifelse(log_Mut_freq_change > 0, log_Mut_freq_change/log_max_change,
                                       -log_Mut_freq_change/log_min_change))

## Identify points for invasion boundary line ----
sweep0_d <- sweep_plot_f %>%
  filter(Mut_freq_inv == "Decrease") %>%
  slice_min(order_by = abs(log_Mut_freq_change2), n = 100)
sweep0_i <- sweep_plot_f %>%
  filter(Mut_freq_inv == "Increase") %>%
  slice_min(order_by = abs(log_Mut_freq_change2), n = 100)
sweep0 <- rbind(sweep0_d, sweep0_i)


# Final plot ----
plot <- ggplot() +
  geom_raster(data = sweep_plot_f,
              mapping = aes(log_gamma_M, psi_M, fill = log_Mut_freq_change2)) +
  scale_fill_gradient2("Relative\nfrequency\nchange",
                       low=p_highD, mid = p_mid, high=p_highI, midpoint = log(1)) +
  geom_hline(yintercept = psi_ref, color = p_axes, linewidth = 1) +
  geom_vline(xintercept = gamma_ref, color = p_axes, linewidth = 1) +
  geom_contour(data = sweep_plot_f,
               mapping = aes(x = log_gamma_M, psi_M, z = log_Mut_freq_change2),
               breaks = 0, color = p_invline, linewidth = inv_width) +
  labs(x = expression("log10(Conjugation Rate)"),
       y = expression("Growth Rate")) +
  scale_x_continuous(expand = c(0.025, 0.025),
                     limits = c(min_gamma, max_gamma)) +
  scale_y_continuous(expand = c(0.01, 0.01),
                     limits = c(min_psi, max_psi)) +
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
  fig_aes +
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
