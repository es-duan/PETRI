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

## Points ----
plot_points <- jsonlite::fromJSON(args$points)
psweep_point_size <- plot_points[["psweep_point_size"]]
sh_Anc <- plot_points[["sh_Anc"]]
sh_Mut <- plot_points[["sh_Mut"]]
sh_SynMut <- plot_points[["sh_SynMut"]]

## Lines ----
plot_lines <- fromJSON(args$lines)
inv_width <- plot_lines[["inv_width"]]

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
  geom_raster(data = filter(p1_DG2, Mut_freq_inv == "Increase"),
              mapping = aes(log_gamma_M, psi_M, fill = log_Mut_freq_change2)) +
  geom_raster(data = filter(p1_DG2, Mut_freq_inv == "Decrease"),
              mapping = aes(log_gamma_M, psi_M),
              fill = p_highD) +
  scale_fill_gradient("Relative\nfrequency\nchange",
                       low = p_mid, high=p_highI) +
  geom_hline(yintercept = psi_ref, color = p_axes, linewidth = 1) +
  geom_vline(xintercept = log10(gamma_ref), color = p_axes, linewidth = 1) +
  geom_point(data = ph_s,
             mapping = aes(log_conj, Growth_rate, shape = Strain),
             color = "black", size = psweep_point_size) +
  scale_shape_manual(values = c("S.pB10" = sh_Anc,
                                "S.pB10-A" = sh_Mut,
                                "S.pB10-B" = sh_SynMut),
                     labels = c("S.pB10" = "Anc",
                                "S.pB10-A" = "Mut",
                                "S.pB10-B" = "Syn Mut"),
                     name = "Genotype") +
  # geom_contour(data = p1_DG2,
  #              mapping = aes(x = log_gamma_M, psi_M, z = log_Mut_freq_change2),
  #              breaks = 0, color = p_invline, linewidth = inv_width) +
  # geom_smooth(data = sweep0_p1,
  #             mapping = aes(log_gamma_M, psi_M),
  #             color = p_invline, method = "lm", se = FALSE, linewidth = inv_width) +
  annotate(geom = "text", x = -14.75, y = 0.1,
           label = "exclusion", color = "white", fontface = "italic",
           size = 4) +
  scale_x_continuous(expand = c(0.01, 0.01)) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  labs(x = expression("log10(Conjugation Rate)"),
       y = expression("Growth Rate")) +
  fig_aes +
  guides(shape = guide_legend(order = 1),
         fill = guide_colourbar(display = "gradient",
                                ticks = FALSE,
                                theme = theme(legend.text = element_text(size = 6)),
                                barwidth = unit(0.15, "in"), barheight = unit(0.6, "in")))

plot1n <- plot1 + 
  theme(legend.position = "none")


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
  geom_raster(data = filter(p2_DGs2, Mut_freq_inv == "Increase"),
              mapping = aes(log_gamma_M, psi_M, fill = log_Mut_freq_change2)) +
  geom_raster(data = filter(p2_DGs2, Mut_freq_inv == "Decrease"),
              mapping = aes(log_gamma_M, psi_M),
              fill = p_highD) +
  scale_fill_gradient("Relative\nfrequency\nchange",
                       low = p_mid, high=p_highI) +
  geom_hline(yintercept = psi_ref, color = p_axes, linewidth = 1) +
  geom_vline(xintercept = log10(gamma_ref), color = p_axes, linewidth = 1) +
  geom_point(data = ph_s,
             mapping = aes(log_conj, Growth_rate, shape = Strain),
             color = "black", size = psweep_point_size) +
  scale_shape_manual(values = c("S.pB10" = sh_Anc,
                                "S.pB10-A" = sh_Mut,
                                "S.pB10-B" = sh_SynMut),
                     labels = c("S.pB10" = "Anc",
                                "S.pB10-A" = "Mut",
                                "S.pB10-B" = "Syn Mut"),
                     name = "Genotype") +
  # geom_smooth(data = sweep0_p2,
  #             mapping = aes(log_gamma_M, psi_M),
  #             color = p_invline, method = "lm", se = FALSE, linewidth = inv_width) +
  annotate(geom = "text", x = -14.75, y = 0.1,
           label = "exclusion", color = "white", fontface = "italic",
           size = 4) +
  scale_x_continuous(expand = c(0.01, 0.01)) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  labs(x = expression("log10(Conjugation Rate)"),
       y = expression("Growth Rate")) +
  fig_aes +
  guides(shape = guide_legend(order = 1),
         fill = guide_colourbar(display = "gradient",
                                ticks = FALSE,
                                theme = theme(legend.text = element_text(size = 6)),
                                barwidth = unit(0.15, "in"), barheight = unit(0.6, "in")))


final_plot <- plot1n + plot2 +
  plot_layout(nrow = 1,
              axes = "collect",
              axis_titles = "collect",
              guides = "collect")

ggsave("figures/figS3_DGsweeps.pdf",
       final_plot, width = 7, height = 3.25, units = "in")
