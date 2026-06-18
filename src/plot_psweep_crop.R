# Crop and plot main psweeps

# Load packages ----
library(tidyverse)
library(scales)
library(argparse)
library(jsonlite)

# Set arguments parser inputs ----
parser <- ArgumentParser()
parser$add_argument("-p","--psweepsetting", type = "character", help = "Specify Treatment and parameter sweep setting")
parser$add_argument("-c","--colors", help = "JSON string of plot colors")
parser$add_argument("-l","--lines", help = "JSON string of plot lines")

# Parse arguments
args <- parser$parse_args()

# Get treatment
ps <- args$psweepsetting
output_folder <- paste("results/parameter_sweeps", ps, sep = "/")

# Read in files ----
sweep_plot <- read_csv(paste0(output_folder, "/" , ps, "_plot.csv"))
setting_list <- readRDS(paste0(output_folder, "/", ps, "_settings.rds"))
ph <- read_csv("input_data/strain_phenotypes.csv") %>%
  mutate(log_Conj = log10(Conjugation_rate)) %>%
  filter(!is.na(Strain))
frequency_limits <- readRDS("figures/panels/fig6_limits.rds")

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

## Lines ----
plot_lines <- fromJSON(args$lines)
inv_width <- plot_lines[["inv_width"]]

## Retrieve ggplot theme ----
source("src/ggplot_theme.R")

# Plot psweep ----

## Set plot limits ----
anc_gamma <- ph$log_Conj[ph$Strain == "S.pB10"]
mut_gamma <- ph$log_Conj[ph$Strain == "S.pB10-A"]
anc_psi <- ph$Growth_rate[ph$Strain == "S.pB10"]
mut_psi <- ph$Growth_rate[ph$Strain == "S.pB10-A"]

mid_gamma <- mean(c(anc_gamma, mut_gamma))
mid_psi <- mean(c(anc_psi, mut_psi))

max_gamma <- mid_gamma + 3.5
min_gamma <- mid_gamma - 3.5

max_psi <- mid_psi + 0.35
min_psi <- mid_psi - 0.35

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


# Plot with geom smooth line ----
i2 <- ggplot() +
  geom_raster(data = sweep_plot_f,
              mapping = aes(log_gamma_M, psi_M, fill = log_Mut_freq_change2)) +
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
  scale_x_continuous(limits = c(min_gamma, max_gamma),
                     expand = c(0.025, 0.025)) +
  scale_y_continuous(limits = c(min_psi, max_psi),
                     expand = c(0.01, 0.01)) +
  labs(x = expression("log10(Conjugation Rate)"),
       y = expression("Growth Rate")) +
  fig_aes 

# Save plot
# ggsave(paste0(output_folder, "/", ps, "_inv_change_plot2_crop.pdf"),
#        i2, height = 2.5, width = 3.5, units = "in")
saveRDS(i2, paste0(output_folder, "/", ps, "_inv_change_plot2_crop.rds"))

### Geom smooth line with no scaling ----
i2n <- ggplot() +
  geom_raster(data = sweep_plot_f,
              mapping = aes(log_gamma_M, psi_M, fill = log_Mut_freq_change)) +
  scale_fill_gradientn("log10\nFrequency\nchange",
                       colors = c(p_highD, p_mid, p_highI),
                       values = scales::rescale(c(frequency_limits[1], log(1), frequency_limits[2])),
                       limits = frequency_limits,
                       breaks = c(-18, -9, 0, 2),
                       oob = scales::squish) +
  geom_hline(yintercept = psi_ref, color = p_axes, linewidth = 1) +
  geom_vline(xintercept = log10(gamma_ref), color = p_axes, linewidth = 1) +
  geom_smooth(data = sweep0,
              mapping = aes(log_gamma_M, psi_M),
              color = p_invline, se = FALSE, linewidth = inv_width) +
  # geom_text(data = axes_label,
  #           mapping = aes(x, y, label = label),
  #           size = 4) +
  scale_x_continuous(limits = c(min_gamma, max_gamma),
                     expand = c(0.025, 0.025)) +
  scale_y_continuous(limits = c(min_psi, max_psi),
                     expand = c(0.01, 0.01)) +
  labs(x = expression("log10(Conjugation Rate)"),
       y = expression("Growth Rate")) +
  fig_aes 

# Save plot
# ggsave(paste0(output_folder, "/", ps, "_inv_change_plot2_ns.pdf"),
#        i2n, height = 2.5, width = 3.5, units = "in")
saveRDS(i2n, paste0(output_folder, "/", ps, "_inv_change_plot2_ns_crop.rds"))
