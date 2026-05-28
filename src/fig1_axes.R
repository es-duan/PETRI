# Draw axes plot

# Load packages ----
library(tidyverse)
library(funkyheatmap)
library(argparse)
library(jsonlite)

# Set arguments parser inputs ----
parser <- ArgumentParser()
parser$add_argument("-t","--treatment", type = "character", help = "Specify Treatment ID")
parser$add_argument("-c","--colors", help = "JSON string of plot colors")

# Parse arguments
args <- parser$parse_args()

# Load global variables ----
## Colors ----
plot_colors <- fromJSON(args$colors)
p_hc <- plot_colors[["p_hc"]]
p_syn <- plot_colors[["p_syn"]]
p_par <- plot_colors[["p_par"]]

## Retrieve ggplot theme ----
source("src/ggplot_theme.R")

# Data frame for rectangles ----
rect_df <- data.frame("Quadrant" = c("Host-centric", "Synergistic", "Parasitic"),
                      "xmin" = c(-10, 0, 0),
                      "xmax" = c(0, 10, 10),
                      "ymin" = c(0, 0, -10),
                      "ymax" = c(10, 10, 0))

rect_df2 <- data.frame("Quadrant" = rep(c("Host-centric", "Synergistic", "Parasitic"),2),
                      "xmin" = c(-10, 0, 0, -1, 0, 0),
                      "xmax" = c(0, 10, 10, 0, 1, 1),
                      "ymin" = c(0, 0, -1, 9, 9, -9),
                      "ymax" = c(1, 1, 0, 10, 10, -10))


# ggplot ----
p1 <- ggplot() +
  geom_rounded_rect(mapping = aes(xmin = -10, xmax = 10, ymin = -10, ymax = 10),
                    color = "gray70", fill = "white", radius = unit(3, "mm")) +
  geom_rounded_rect(mapping = aes(xmin = -9.5, xmax = -1, ymin = -9.5, ymax = -5.5),
                    color = "gray70", fill = "white", radius = unit(3, "mm")) +
  geom_rect(data = rect_df2,
                    mapping = aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = Quadrant)) +
  geom_rounded_rect(data = rect_df,
                    mapping = aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = Quadrant),
                    radius = unit(3, "mm")) +
  geom_segment(mapping = aes(x = -10.2, xend = 10.2, y = 0, yend = 0),
               color = "black",
               arrow = arrow(length = unit(0.2, "cm"),
                             type = "closed",
                             ends = "both"),
               linewidth = 1,
               lineend = "round") +
  geom_segment(mapping = aes(x = -0, xend = 0, y = -10.2, yend = 10.2),
               color = "black",
               arrow = arrow(length = unit(0.2, "cm"),
                             type = "closed",
                             ends = "both"),
               linewidth = 1,
               lineend = "round") +
  geom_point(mapping = aes(x = 0, y = 0),
             color = "black", size = 3) +
  scale_x_continuous(limits = c(-10.5,10.5),
                     name = expression(paste(Delta, "Conjugation (HGT)")),
                     expand = c(0.01, 0.01)) +
  scale_y_continuous(limits = c(-10.5,10.5),
                     name = expression(paste(Delta, "Host Growth (VGT)")),
                     expand = c(0.01, 0.01)) +
  scale_fill_manual(values = c("Host-centric" = p_hc,
                               "Synergistic" = p_syn,
                               "Parasitic" = p_par),
                    name = NULL) +
  fig_aes +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        legend.position = c(0.25, 0.15),
        legend.key.size = unit(0.1, "in"),
        legend.background = element_rect(fill = "transparent", color = NA))

ggsave("figures/panels/fig1_axes.pdf",
       p1, width = 3, height = 3, units = "in")
saveRDS(p1, "figures/panels/fig1_axes.rds")

