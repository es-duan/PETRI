# Fig S4: replot/clean up Dim copA invasion

# Load packages ----
library(tidyverse)
library(patchwork)
library(argparse)
library(jsonlite)

# Set arguments parser inputs ----
parser <- ArgumentParser()
parser$add_argument("-c","--colors", help = "JSON string of plot colors")
parser$add_argument("-l","--lines", help = "JSON string of plot lines")

# Parse arguments
args <- parser$parse_args()

# Read in files ----
sim_dens <- read_csv("results/case_study_sims/Dim90_copA/Dim90_copA_density_plot_df.csv")
sim_freq <- read_csv("results/case_study_sims/Dim90_copA/Dim90_copA_frequency_plot_df.csv")
phases_rect <- read_csv("results/case_study_sims/Dim90_copA/Dim90_copA_phases_plot_df.csv")

# Load global variables ----
## Colors ----
plot_colors <- fromJSON(args$colors)
p_Anc <- plot_colors[["p_Anc"]]
p_Mut <- plot_colors[["p_Mut"]]
p_F <- plot_colors[["p_F"]]
p_growth <- plot_colors[["p_growth"]]
p_conj <- plot_colors[["p_conj"]]
p_tselect <- plot_colors[["p_tselect"]]
p_imm <- plot_colors[["p_imm"]]

## Lines ----
plot_lines <- fromJSON(args$lines)
rifR_l <- plot_lines[["rifR_l"]]
nalR_l <- plot_lines[["nalR_l"]]
plot_lw <- plot_lines[["plot_lw"]]

## Retrieve ggplot theme ----
source("src/ggplot_theme.R")

# Edit density plot ----
pA <- ggplot() + 
  geom_rect(data = phases_rect,
            mapping = aes(xmin = T_start, xmax = T_end, ymin = 1, ymax = Inf, fill = Phase)) +
  scale_fill_manual(values = c("Growth" = p_growth,
                               "Conjugation" = p_conj,
                               "Transconjugant selection" = p_tselect,
                               "Immigration" = p_imm),
                    labels = c("Immigration" = "90% immigration")) +
  geom_line(data = sim_dens,
            mapping = aes(x = Time, y = Density,
                          color = Genotype_plot, linetype = Host,
                          group = interaction(Cycle, Genotype_plot, Host)),
            linewidth = plot_lw) +
  scale_color_manual(values = c("Ancestor" = p_Anc,
                                "Mutant" = p_Mut,
                                "Plasmid-free" = p_F),
                     breaks = c("Ancestor", "Mutant", "Plasmid-free"),
                     labels = c("R1", "copA", "Plasmid-free")) +
  scale_linetype_manual(values = c("rifR" = rifR_l,
                                   "nalR" = nalR_l)) +
  #geom_hline(yintercept = 1.5, color = "white") +
  labs(x = "Time (hr)",
       y = "Density (CFU/mL)",
       color = "Genotype",
       fill = "Phase",
       linetype = "Host") +
  scale_y_continuous(trans = "log10",
                     breaks=c(1e0,1e2,1e4,1e6,1e8),
                     labels=sapply(c(0,2,4,6,8),function(i){parse(text = sprintf("10^%d",i))}),
                     expand = c(0.01, 0.01),
                     limits = c(0.9, 1.5e9)) +
  scale_x_continuous(expand = c(0.001, 0.001)) +
  guides(color = guide_legend(order=1, override.aes = list(shape = NA)),
         fill = guide_legend(order=2),
         linetype = "none") +
  fig_aes +
  theme(legend.position = "top",
        legend.title = element_blank())

# Edit frequency plot ----
sim_freq2 <- sim_freq %>%
  mutate(Frequency = ifelse(is.na(Frequency),0,Frequency))

f_limits <- c(1e-4, 1)
f_breaks <- c(1e-4, 1e-2, 1)
f_labels <- c(-4, -2)

pB <- ggplot() + 
  geom_rect(data = phases_rect,
            mapping = aes(xmin = T_start, xmax = T_end, ymin = 0, ymax = Inf, fill = Phase)) +
  scale_fill_manual(values = c("Growth" = p_growth,
                               "Conjugation" = p_conj,
                               "Transconjugant selection" = p_tselect,
                               "Immigration" = p_imm)) +
  geom_line(data = sim_freq2,
            mapping = aes(Time, Frequency, color = Genotype),
            linewidth = plot_lw) +
  scale_color_manual(values = c("Ancestor" = p_Anc,
                                "Mutant" = p_Mut)) +
  scale_y_continuous(limits = f_limits,
                     breaks = f_breaks,
                     labels = c(sapply(f_labels,function(i){parse(text = sprintf("10^%d",i))}),1),
                     expand = c(0.01, 0.01),
                     trans = "log10") +
  scale_x_continuous(expand = c(0.001, 0.001)) +
  labs(x = "Time (hr)",
       y = "Frequency",
       color = "Genotype",
       fill = "Phase") +
  fig_aes +
  theme(legend.position = "none")

# Combine plots ----
final_plot <- pA + pB +
  plot_layout(nrow = 2,
              axes = "collect",
              axis_titles = "collect") +
  plot_annotation(tag_levels = "A")

# Save plot
ggsave("figures/figS4_copA_sim.pdf",
       final_plot, width = 7, height = 5, units = "in")

