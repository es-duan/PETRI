# Supplementary parameter sweep plots

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
p_Exc <- plot_colors[["p_Exc"]]
p_Dis <- plot_colors[["p_Dis"]]
p_lowI <- plot_colors[["p_lowI"]]
p_highI <- plot_colors[["p_highI"]]
p_lowD <- plot_colors[["p_lowD"]]
p_highD <- plot_colors[["p_highD"]]
p_axes <- plot_colors[["p_axes"]]
p_mid <- plot_colors[["p_mid"]]
p_invline <- plot_colors[["p_invline"]]

inv_width <- 0.7

## Points ----
plot_points <- jsonlite::fromJSON(args$points)
psweep_point_size <- plot_points[["psweep_point_size"]]
sh_Anc <- plot_points[["sh_Anc"]]
sh_Mut <- plot_points[["sh_Mut"]]
sh_R1 <- plot_points[["sh_R1"]]
sh_copA <- plot_points[["sh_copA"]]
sh_finO <- plot_points[["sh_finO"]]

## Retrieve ggplot theme ----
source("src/ggplot_theme.R")

## Read in files ----
p1_DG <- read_csv("results/parameter_sweeps/pDG_S.pB10/pDG_S.pB10_plot.csv")
p2_DGs <- read_csv("results/parameter_sweeps/pDGs_S.pB10/pDGs_S.pB10_plot.csv")

ph <- read_csv("input_data/strain_phenotypes.csv")

# Select relevant strains to plot ----
ph_s <- ph %>%
  filter(str_detect(Strain, "S.")) %>%
  #filter(Strain != ps_s) %>%
  mutate(log_conj = log10(Conjugation_rate))

gamma_ref <- ph_s$Conjugation_rate[ph_s$Strain == "S.pB10"]
psi_ref <- ph_s$Growth_rate[ph_s$Strain == "S.pB10"]

# Geom smooth line ----
## P1 ----
log_max_change1 <- max(p1_DG$log_Mut_freq_change)
log_min_change1 <- min(p1_DG$log_Mut_freq_change)

p1_DG2 <- p1_DG %>%
  mutate(log_Mut_freq_change2 = case_when(Mut_freq == 0 ~ -1,
                                          log_Mut_freq_change > 0 ~ log_Mut_freq_change/log_max_change1,
                                          log_Mut_freq_change <= 0 ~ -log_Mut_freq_change/log_min_change1))

### Identify points for invasion boundary line ----
sweep0_p1 <- p1_DG2 %>% 
  filter(abs(log_Mut_freq_change2) < 0.1)

plot1 <- ggplot() +
  geom_raster(data = p1_DG2,
              mapping = aes(log_gamma_M, psi_M, fill = log_Mut_freq_change2)) +
  scale_fill_gradient2("Relative\nfrequency\nchange",
                       low=p_highD, mid = p_mid, high=p_highI, midpoint = log(1)) +
  geom_hline(yintercept = psi_ref, color = p_axes, linewidth = 1) +
  geom_vline(xintercept = log10(gamma_ref), color = p_axes, linewidth = 1) +
  geom_point(data = ph_s,
             mapping = aes(log_conj, Growth_rate, shape = Strain),
             color = "black", size = psweep_point_size) +
  scale_shape_manual(values = c("S.pB10" = sh_Anc,
                                "S.pB10-A" = sh_Mut,
                                "S.pB10-B" = sh_finO)) +
  # geom_smooth(data = sweep0_p1,
  #             mapping = aes(log_gamma_M, psi_M),
  #             color = p_invline, method = "lm", se = FALSE, linewidth = inv_width) +
  scale_x_continuous(expand = c(0.01, 0.01)) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  labs(x = expression("log10(Conjugation Rate)"),
       y = expression("Growth Rate")) +
  fig_aes 


## P2 ----
log_max_change2 <- max(p2_DGs$log_Mut_freq_change)
log_min_change2 <- min(p2_DGs$log_Mut_freq_change)

p2_DGs2 <- p2_DGs %>%
  mutate(log_Mut_freq_change2 = case_when(Mut_freq == 0 ~ -1,
                                          log_Mut_freq_change > 0 ~ log_Mut_freq_change/log_max_change2,
                                          log_Mut_freq_change <= 0 ~ -log_Mut_freq_change/log_min_change2))

### Identify points for invasion boundary line ----
sweep0_p2 <- p2_DGs2 %>% 
  filter(abs(log_Mut_freq_change2) < 0.1)

plot2 <- ggplot() +
  geom_raster(data = p2_DGs2,
              mapping = aes(log_gamma_M, psi_M, fill = log_Mut_freq_change2)) +
  scale_fill_gradient2("Relative\nfrequency\nchange",
                       low=p_highD, mid = p_mid, high=p_highI, midpoint = log(1)) +
  geom_hline(yintercept = psi_ref, color = p_axes, linewidth = 1) +
  geom_vline(xintercept = log10(gamma_ref), color = p_axes, linewidth = 1) +
  geom_point(data = ph_s,
             mapping = aes(log_conj, Growth_rate, shape = Strain),
             color = "black", size = psweep_point_size) +
  scale_shape_manual(values = c("S.pB10" = sh_Anc,
                                "S.pB10-A" = sh_Mut,
                                "S.pB10-B" = sh_finO)) +
  # geom_smooth(data = sweep0_p2,
  #             mapping = aes(log_gamma_M, psi_M),
  #             color = p_invline, method = "lm", se = FALSE, linewidth = inv_width) +
  scale_x_continuous(expand = c(0.01, 0.01)) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  labs(x = expression("log10(Conjugation Rate)"),
       y = expression("Growth Rate")) +
  fig_aes 


final_plot <- plot1 + plot2 +
  plot_layout(nrow = 1,
              axes = "collect",
              axis_titles = "collect",
              guides = "collect") +
  plot_annotation(tag_levels = "a")

ggsave("figures/fig_s1_DGsweeps.pdf",
       final_plot, width = 10.5, height = 5, units = "in")
