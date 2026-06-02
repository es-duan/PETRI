# Plot parameter sweep data

# Load packages ----
library(tidyverse)
library(argparse)
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
setting_list <- readRDS(paste0(output_folder, "/", ps, "_settings.rds"))
ph <- read_csv("input_data/strain_phenotypes.csv")

# Relevant parameters
gamma_ref <- as.numeric(setting_list$Ref_gamma)
psi_ref <- as.numeric(setting_list$Ref_psi)

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

## Retrieve ggplot theme ----
source("src/ggplot_theme.R")

# Plot psweep ----

## Set plot limits ----
gamma_M_max = max(sweep_plot$log_gamma_M)
gamma_M_min = min(sweep_plot$log_gamma_M)
psi_M_max = max(sweep_plot$psi_M)
psi_M_min = min(sweep_plot$psi_M)

## Select axes labels to plot ----
gamma_all = seq(round(gamma_M_min, 0), round(gamma_M_max, 0), 1)
gamma_range = gamma_all[seq(2,length(gamma_all), 2)]

psi_all = seq(round(psi_M_min, 1), round(psi_M_max, 1), 0.1)
psi_range = psi_all[seq(2, length(psi_all), 2)]

# Generate dataset for axes labels
axes1 <- data.frame("gamma" = c(gamma_range, rep(log10(gamma_ref), length(psi_range))),
                      "psi" = c(rep(psi_ref, length(gamma_range)), psi_range)) %>%
  mutate(x = ifelse(gamma == log10(gamma_ref), gamma - 0.25, gamma),
         y = ifelse(psi == psi_ref, psi - 0.015, psi)) %>%
  mutate(label = case_when(gamma == log10(gamma_ref) ~ as.character(psi),
                           psi == psi_ref ~ as.character(gamma)))

axes2 <- data.frame("gamma" = log10(gamma_ref),
                    "psi" = psi_ref) %>%
  mutate(x = gamma + 0.8,
         y = psi + 0.03) %>%
  mutate(label = paste0("(",as.character(round(gamma, 1)),",",
                        as.character(round(psi,2)),")"))

axes_label <- rbind(axes1,axes2) %>%
  mutate(label_f = case_when(label == axes2$label[1] ~ label,
                             psi == psi_ref & gamma != gamma_ref ~ as.character(round(log10(10^gamma/gamma_ref), 2)),
                             gamma == log10(gamma_ref) & psi != psi_ref ~ as.character(round(psi - psi_ref, 2))))

## Transformations for color ----
log_max_change <- max(sweep_plot$log_Mut_freq_change)
log_min_change <- min(sweep_plot$log_Mut_freq_change)

sweep_plot2 <- sweep_plot %>%
  mutate(log_Mut_freq_change2 = ifelse(log_Mut_freq_change > 0, log_Mut_freq_change/log_max_change,
                                       -log_Mut_freq_change/log_min_change))

## Identify points for invasion boundary line ----
sweep0_d <- sweep_plot2 %>%
  filter(Mut_freq_inv == "Decrease") %>%
  slice_min(order_by = abs(log_Mut_freq_change2), n = 100)
sweep0_i <- sweep_plot2 %>%
  filter(Mut_freq_inv == "Increase") %>%
  slice_min(order_by = abs(log_Mut_freq_change2), n = 100)
sweep0 <- rbind(sweep0_d, sweep0_i)

# Plot by relative frequency change ----
## Contour line ----
i1 <- ggplot() +
  geom_raster(data = sweep_plot2,
              mapping = aes(log_gamma_M, psi_M, fill = log_Mut_freq_change2)) +
  scale_fill_gradient2("Relative\nfrequency\nchange",
                       low=p_highD, mid = p_mid, high=p_highI, midpoint = log(1)) +
  geom_hline(yintercept = psi_ref, color = p_axes, linewidth = 1) +
  geom_vline(xintercept = log10(gamma_ref), color = p_axes, linewidth = 1) +
  geom_contour(data = sweep_plot2,
               mapping = aes(x = log_gamma_M, psi_M, z = log_Mut_freq_change2),
               breaks = 0, color = p_invline, linewidth = inv_width) +
  # geom_text(data = axes_label,
  #           mapping = aes(x, y, label = label),
  #           size = 4) +
  scale_x_continuous(expand = c(0.01, 0.01)) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  labs(x = expression("log10(Conjugation Rate)"),
       y = expression("Growth Rate")) +
  fig_aes 

# Save plot
ggsave(paste0(output_folder, "/", ps, "_inv_change_plot.pdf"),
       i1, height = 2.5, width = 3.5, units = "in")
saveRDS(i1, paste0(output_folder, "/", ps, "_inv_change_plot.rds"))

# Geom smooth line ----
i2 <- ggplot() +
  geom_raster(data = sweep_plot2,
              mapping = aes(log_gamma_M, psi_M, fill = log_Mut_freq_change2)) +
  # geom_contour(data = sweep_plot2,
  #              mapping = aes(x = log_gamma_M, psi_M, z = log_Mut_freq_change2),
  #              breaks = 0, color = p_invline) +
  scale_fill_gradient2("Relative\nfrequency\nchange",
                       low=p_highD, mid = p_mid, high=p_highI, midpoint = log(1)) +
  geom_hline(yintercept = psi_ref, color = p_axes, linewidth = 1) +
  geom_vline(xintercept = log10(gamma_ref), color = p_axes, linewidth = 1) +
  geom_smooth(data = sweep0,
              mapping = aes(log_gamma_M, psi_M),
              color = p_invline, se = FALSE, linewidth = inv_width) +
  # geom_text(data = axes_label,
  #           mapping = aes(x, y, label = label),
  #           size = 4) +
  scale_x_continuous(expand = c(0.01, 0.01)) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  labs(x = expression("log10(Conjugation Rate)"),
       y = expression("Growth Rate")) +
  fig_aes 

# Save plot
ggsave(paste0(output_folder, "/", ps, "_inv_change_plot2.pdf"),
       i2, height = 2.5, width = 3.5, units = "in")
saveRDS(i2, paste0(output_folder, "/", ps, "_inv_change_plot2.rds"))

## Binary plot (increase or decrease) ----
i2 <- ggplot() +
  geom_tile(data = sweep_plot,
            aes(log_gamma_M, psi_M, fill = Mut_freq_inv),
            color = NA) +
  geom_hline(yintercept = psi_ref, color = "black", linewidth = 1) + 
  geom_vline(xintercept = log10(gamma_ref), color = "black", linewidth = 1) +
  geom_text(data = axes_label,
            mapping = aes(x, y, label = label),
            size = 4) +
  scale_fill_manual(values = c("Increase" = p_Dis,
                               "Decrease" = p_Exc,
                               "NA" = "gray")) +
  scale_x_continuous(expand = c(0.01, 0.01),
                     labels = ~ ifelse(.x == 0, "", .x)) +
  scale_y_continuous(expand = c(0.01, 0.01),
                     labels = ~ ifelse(.x == 0, "", .x)) +
  labs(x = expression("log10(Conjugation Rate)"),
       y = expression("Growth Rate"),
       fill = "Invasion") +
  fig_aes +
  theme(axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# Save plot
ggsave(paste0(output_folder, "/", ps, "_change_plot.pdf"),
       i2, height = 2.5, width = 3.5, units = "in")
