# Process parameter sweep data

# Load packages ----
library(tidyverse)
library(argparse)
library(ggh4x)
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
sweep_plot <- read_csv(paste0(output_folder, "/" , ps, "_plot.csv"))

# Load global variables ----
## Colors ----
plot_colors <- jsonlite::fromJSON(args$colors)
p_Exc <- plot_colors[["p_Exc"]]
p_Dis <- plot_colors[["p_Dis"]]
p_lowI <- plot_colors[["p_lowI"]]
p_highI <- plot_colors[["p_highI"]]
p_lowD <- plot_colors[["p_lowD"]]
p_highD <- plot_colors[["p_highD"]]

# Non-snakemake
# p_Exc <- "gray95"
# p_Dis <- "#140433"
# p_lowD <- "#BABAFF"
# p_highD <- "#0000AB"
# p_lowI <- "#FFC2C2"
# p_highI <- "#C20000"

## Retrieve ggplot theme ----
source("src/ggplot_theme.R")

# Plot psweep ----

# Set plot limits
gamma_fold_max = max(sweep_plot$gamma_fold)
gamma_fold_min = min(sweep_plot$gamma_fold)
psi_change_max = max(sweep_plot$psi_change)
psi_change_min = min(sweep_plot$psi_change)

# Generate dataset for axes labels
axes_label <- sweep_plot %>%
  filter(gamma_fold %% 1 == 0) %>%
  filter(psi_change %% 0.1 == 0) %>%
  filter(gamma_fold == 0 | psi_change == 0) %>%
  filter(!(gamma_fold == 0 & psi_change == 0)) %>%
  select(gamma_fold, psi_change) %>%
  mutate(x = ifelse(gamma_fold == 0, gamma_fold - 0.25, gamma_fold),
         y = ifelse(psi_change == 0, psi_change - 0.015, psi_change)) %>%
  mutate(label = ifelse(gamma_fold == 0, psi_change, gamma_fold))

## Plot by rate of change ----
i1 <- ggplot() +
  geom_tile(data = filter(sweep_plot, Mut_freq_inv == "Increase"),
            mapping = aes(gamma_fold,psi_change, fill = Mut_freq_change)) +
  scale_fill_gradient("Frequency\nincrease",
                      low=p_lowI, high=p_highI) +
  new_scale_fill() +
  geom_tile(data = filter(sweep_plot, Mut_freq_inv == "Decrease"),
            mapping = aes(gamma_fold,psi_change, fill = Mut_freq_change)) +
  scale_fill_gradient("Frequency\ndecrease",
                      low=p_highD, high=p_lowD) +
  geom_hline(yintercept = 0, color = "black", linewidth = 1.5) + 
  geom_vline(xintercept = 0, color = "black", linewidth = 1.5) +
  geom_text(data = axes_label,
            mapping = aes(x, y, label = label),
            size = 6) +
  scale_x_continuous(expand = c(0.01, 0.01),
                     labels = ~ ifelse(.x == 0, "", .x)) +
  scale_y_continuous(expand = c(0.01, 0.01),
                     labels = ~ ifelse(.x == 0, "", .x)) +
  coord_axes_inside(labels_inside = TRUE,
                    xlim=c(gamma_fold_min,gamma_fold_max),
                    ylim=c(psi_change_min,psi_change_max)) +
  labs(x = expression(Delta*"Conjugation Rate"),
       y = expression(Delta*"Host Fitness")) +
  axes_aes

# Save plot
ggsave(paste0(output_folder, "/", ps, "_inv_change_plot.pdf"),
       i1, height = 5, width = 6.75, units = "in")

## Binary plot (increase or decrease) ----
i2 <- ggplot() +
  geom_tile(data = sweep_plot,
            aes(gamma_fold, psi_change, fill = Mut_freq_inv)) +
  geom_hline(yintercept = 0, color = "black", linewidth = 1.5) + 
  geom_vline(xintercept = 0, color = "black", linewidth = 1.5) +
  geom_text(data = axes_label,
            mapping = aes(x, y, label = label),
            size = 6) +
  scale_fill_manual(values = c("Increase" = p_Dis,
                               "Decrease" = p_Exc,
                               "NA" = "gray")) +
  scale_x_continuous(expand = c(0.01, 0.01),
                     labels = ~ ifelse(.x == 0, "", .x)) +
  scale_y_continuous(expand = c(0.01, 0.01),
                     labels = ~ ifelse(.x == 0, "", .x)) +
  coord_axes_inside(labels_inside = TRUE,
                    xlim=c(gamma_fold_min,gamma_fold_max),
                    ylim=c(psi_change_min,psi_change_max)) +
  labs(x = expression(Delta*"Conjugation Rate"),
       y = expression(Delta*"Host Fitness"),
       fill = "Invasion") +
  axes_aes

# Save plot
ggsave(paste0(output_folder, "/", ps, "_change_plot.pdf"),
       i2, height = 5, width = 6.75, units = "in")


