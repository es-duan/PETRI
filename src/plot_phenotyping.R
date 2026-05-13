# Plot phenotyping data

# Load packages ----
library(tidyverse)
library(argparse)
library(jsonlite)

# Set arguments parser inputs ----
parser <- ArgumentParser()
parser$add_argument("-c","--colors", help = "JSON string of plot colors")
parser$add_argument("-p","--points", help = "JSON string of point aesthetics")

# Parse arguments
args <- parser$parse_args()

# Load global variables ----
## Colors ----
plot_colors <- fromJSON(args$colors)
p_Anc <- plot_colors[["p_Anc"]]
p_Mut <- plot_colors[["p_Mut"]]
p_F <- plot_colors[["p_F"]]
p_hc <- plot_colors[["p_hc"]]
p_syn <- plot_colors[["p_syn"]]
p_par <- plot_colors[["p_par"]]

## Points ----
plot_points <- jsonlite::fromJSON(args$points)
point_size <- plot_points[["ph_point_size"]]
sh_Anc <- plot_points[["sh_Anc"]]
sh_Mut <- plot_points[["sh_Mut"]]

## Retrieve ggplot theme ----
source("src/ggplot_theme.R")

# Load in data ----
phenotyping <- read_csv("results/phenotyping/phenotyping_av.csv") %>%
  filter(Genotype != "F")

## Set plotting parameters
# Plot centers
c_conj <- 10^mean(log10(phenotyping$Conjugation_rate_mean))
c_growth <- mean(phenotyping$Growth_rate_mean)

# Plot dimensions
conj_step <- 10^2.5
growth_step <- 0.3

lim_conj <- c(c_conj/conj_step, c_conj*conj_step)
lim_growth <- c(c_growth - growth_step, c_growth + growth_step)

# Plot growth and conjugation rates ----
p1 <- ggplot() +
  geom_errorbar(data = phenotyping,
                mapping= aes(x = Conjugation_rate_mean,
                             ymin = Growth_rate_mean - Growth_rate_se,
                             ymax = Growth_rate_mean + Growth_rate_se,
                             color = Genotype),
                width = 0) +
  geom_errorbarh(data = phenotyping,
                 mapping = aes(y = Growth_rate_mean,
                               xmin = Conjugation_rate_mean - Conjugation_rate_se,
                               xmax = Conjugation_rate_mean + Conjugation_rate_se,
                               color = Genotype),
                 height = 0) +
  geom_point(data = phenotyping,
             mapping = aes(Conjugation_rate_mean, Growth_rate_mean, color = Genotype,
                           shape = Genotype),
             size = point_size) +
  scale_x_continuous(trans = "log10", 
                     name = "Conjugation rate",
                     limits = lim_conj,
                     breaks = c(1e-15,1e-13,1e-11),
                     labels = sapply(c(-15,-13,-11),function(i){parse(text = sprintf("10^%d",i))})) +
  scale_y_continuous(name = "Growth Rate",
                     limits = lim_growth) + 
  scale_color_manual(values = c("Anc" = p_Anc,
                                "Mut" = p_Mut)) +
  scale_shape_manual(values = c("Anc" = sh_Anc,
                                "Mut" = sh_Mut)) +
  fig_aes +
  theme(legend.background = element_rect(fill = "white", color = "gray70"),
        legend.position = c(0.2, 0.15))

ggsave("results/phenotyping/pB10_phenotyping.pdf",
       p1, width = 7.5, height = 7, units = "in")
saveRDS(p1, "results/phenotyping/pB10_phenotyping.rds")


# Plot anc and mut as reference strains only ----
## Mut invasion ----
a_x <- phenotyping$Conjugation_rate_mean[phenotyping$Genotype == "Anc"]
a_y <-phenotyping$Growth_rate_mean[phenotyping$Genotype == "Anc"]

c1 <- ggplot() +
  geom_rect(mapping = aes(xmin = lim_conj[1], xmax = a_x, ymin = lim_growth[1], ymax = a_y),
            fill = "white") +
  geom_rect(mapping = aes(xmin = lim_conj[1], xmax = a_x, ymin = a_y, ymax = lim_growth[2]),
            fill = p_hc) +
  geom_rect(mapping = aes(xmin = a_x, xmax = lim_conj[2], ymin = a_y, ymax = lim_growth[2]),
            fill = p_syn) +
  geom_rect(mapping = aes(xmin = a_x, xmax = lim_conj[2], ymin = lim_growth[1], ymax = a_y),
            fill = p_par) +
  geom_segment(mapping = aes(x = lim_conj[1], xend = lim_conj[2], y = a_y, yend = a_y),
               color = "black") +
  geom_segment(mapping = aes(x = a_x, xend = a_x, y = lim_growth[1], yend = lim_growth[2]),
               color = "black") +
  geom_point(data = phenotyping,
             mapping = aes(Conjugation_rate_mean, Growth_rate_mean, color = Genotype,
                           shape = Genotype),
             size = point_size) +
  scale_x_continuous(trans = "log10", 
                     name = "Conjugation rate",
                     limits = lim_conj,
                     expand = c(0.001, 0.001)) +
  scale_y_continuous(name = "Growth Rate",
                     limits = lim_growth,
                     expand = c(0.001, 0.001)) + 
  scale_color_manual(values = c("Anc" = p_Anc,
                                "Mut" = p_Mut)) +
  scale_shape_manual(values = c("Anc" = sh_Anc,
                                "Mut" = sh_Mut)) +
  theme_void() +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))

ggsave("results/phenotyping/pB10_mut_inv.pdf",
       c1, width = 4, height = 4, units = "in")
saveRDS(c1, "results/phenotyping/pB10_mut_inv.rds")

## Anc invasion ----
m_x <- phenotyping$Conjugation_rate_mean[phenotyping$Genotype == "Mut"]
m_y <-phenotyping$Growth_rate_mean[phenotyping$Genotype == "Mut"]

c2 <- ggplot() +
  geom_rect(mapping = aes(xmin = lim_conj[1], xmax = m_x, ymin = lim_growth[1], ymax = m_y),
            fill = "white") +
  geom_rect(mapping = aes(xmin = lim_conj[1], xmax = m_x, ymin = m_y, ymax = lim_growth[2]),
            fill = p_hc) +
  geom_rect(mapping = aes(xmin = m_x, xmax = lim_conj[2], ymin = m_y, ymax = lim_growth[2]),
            fill = p_syn) +
  geom_rect(mapping = aes(xmin = m_x, xmax = lim_conj[2], ymin = lim_growth[1], ymax = m_y),
            fill = p_par) +
  geom_segment(mapping = aes(x = lim_conj[1], xend = lim_conj[2], y = m_y, yend = m_y),
               color = "black") +
  geom_segment(mapping = aes(x = m_x, xend = m_x, y = lim_growth[1], yend = lim_growth[2]),
               color = "black") +
  geom_point(data = phenotyping,
             mapping = aes(Conjugation_rate_mean, Growth_rate_mean, color = Genotype,
                           shape = Genotype),
             size = point_size) +
  scale_x_continuous(trans = "log10", 
                     name = "Conjugation rate",
                     limits = lim_conj,
                     expand = c(0.001, 0.001)) +
  scale_y_continuous(name = "Growth Rate",
                     limits = lim_growth,
                     expand = c(0.001, 0.001)) + 
  scale_color_manual(values = c("Anc" = p_Anc,
                                "Mut" = p_Mut)) +
  scale_shape_manual(values = c("Anc" = sh_Anc,
                                "Mut" = sh_Mut)) +
  theme_void() +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))

ggsave("results/phenotyping/pB10_anc_inv.pdf",
       c2, width = 4, height = 4, units = "in")
saveRDS(c2, "results/phenotyping/pB10_anc_inv.rds")
