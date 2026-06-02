# Figure 2: analytic model results

# Load packages ----
library(tidyverse)
library(patchwork)
library(argparse)
library(jsonlite)

# Set arguments parser inputs ----
parser <- ArgumentParser()
parser$add_argument("-c","--colors", help = "JSON string of plot colors")

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

# Blend colors for lighter shade
color_blend_D <- colorRampPalette(c(p_mid, p_highD))
seq_D <- color_blend_D(5)
p_lhighD <- seq_D[3]

color_blend_I <- colorRampPalette(c(p_mid, p_highI))
seq_I <- color_blend_I(5)
p_lhighI <- seq_I[4]

## Retrieve ggplot theme ----
source("src/ggplot_theme.R")

# VGT favoring plot ----

## VGT favoring parameters ----
psiA_v <- 1
sigma_v <- 0.05
delta <- 0.2
gammaA_v <- 0.01
alphaN_v <- 0.01

n_eq_v <- (psiA_v*sigma_v)/(psiA_v - delta + gammaA_v + alphaN_v)
slope_v <- -(n_eq_v/(1 - sigma_v))

## HGT favoring parameters ----
psiA_h <- 1
sigma_h <- 0.45
delta <- 0.5
gammaA_h <- 0.01
alphaN_h <- 0.001

n_eq_h <- (psiA_h*sigma_h)/(psiA_h - delta + gammaA_h + alphaN_h)
slope_h <- -(n_eq_h/(1 - sigma_h))

# Create data frame for plotting ----
df <- data.frame("X" = seq(-15, 15, 0.5)) %>%
  mutate(Y_v = slope_v*X,
         Y_h = slope_h*X)

pv <- ggplot(data = df,
       mapping = aes(X, Y_v)) +
  geom_ribbon(mapping = aes(ymin = -Inf, ymax = Y_v),
              fill = p_lhighD) +
  geom_ribbon(mapping = aes(ymin = Y_v, ymax = Inf),
              fill = p_lhighI) +
  geom_hline(yintercept = 0, color = "white", linewidth = 1) +
  geom_vline(xintercept = 0, color = "white", linewidth = 1) +
  geom_abline(intercept = 0, slope = slope_v, color = p_invline,
              linewidth = inv_width) +
  annotate(geom = "text", x = -7, y = 9,
           label = "invasion", color = "white", fontface = "italic",
           size = 4) +
  annotate(geom = "text", x = -6.7, y = -9,
           label = "exclusion", color = "white", fontface = "italic",
           size = 4) +
  scale_y_continuous(name = expression(paste(Delta, "Growth Rate"))) +
  scale_x_continuous(name = expression(paste(Delta, "Conjugation Rate"))) +
  coord_cartesian(
    xlim = c(-10, 10), 
    ylim = c(-10, 10), 
    expand = FALSE) +
  fig_aes +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

ph <- ggplot(data = df,
       mapping = aes(X, Y_h)) +
  geom_ribbon(mapping = aes(ymin = -Inf, ymax = Y_h),
              fill = p_lhighD) +
  geom_ribbon(mapping = aes(ymin = Y_h, ymax = Inf),
              fill = p_lhighI) +
  geom_hline(yintercept = 0, color = "white", linewidth = 1) +
  geom_vline(xintercept = 0, color = "white", linewidth = 1) +
  geom_abline(intercept = 0, slope = slope_h, color = p_invline,
              linewidth = inv_width) +
  scale_y_continuous(name = expression(paste(Delta, "Growth Rate"))) +
  scale_x_continuous(name = expression(paste(Delta, "Conjugation Rate"))) +
  coord_cartesian(
    xlim = c(-10, 10), 
    ylim = c(-10, 10), 
    expand = FALSE) +
  fig_aes +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

final_plot <- pv + ph +
  plot_layout(nrow = 1) +
  plot_annotation(tag_levels = list(c("B","C")))

ggsave("figures/panels/fig2bc_criterion.pdf",
       final_plot, width = 5, height = 2.75, units = "in")

