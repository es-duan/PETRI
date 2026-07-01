# Generate experimental validation figure

# Load packages ----
library(tidyverse)
library(patchwork)

# Read in figures ----
# HFC_sim <- readRDS("results/case_study_sims/HFC_S.pB10/HFC_S.pB10_frequency_plot.rds")
# LFC_sim <- readRDS("results/case_study_sims/LFC_S.pB10-A/LFC_S.pB10-A_frequency_plot.rds")
# HFC_v <- readRDS("results/experimental_validation/HFC/HFC_frequency_plot.rds")
# LFC_v <- readRDS("results/experimental_validation/LFC/LFC_frequency_plot.rds")

# Condensed phases ----
## Read in plots ----
HFC <- readRDS("figures/panels/fig5d_HFC.rds")
LFC <- readRDS("figures/panels/fig5b_LFC.rds")

## Set plotting parameters ----
phenotyping <- read_csv("results/phenotyping/phenotyping_av.csv") %>%
  filter(Genotype != "F")

# Plot centers
c_conj <- 10^mean(log10(phenotyping$Conjugation_rate_mean))
c_growth <- mean(phenotyping$Growth_rate_mean)

# Plot dimensions
conj_step <- 10^2.5
growth_step <- 0.3

lim_conj <- c(c_conj/conj_step, c_conj*conj_step)
lim_growth <- c(c_growth - growth_step, c_growth + growth_step)

a_x <- phenotyping$Conjugation_rate_mean[phenotyping$Genotype == "Anc"]
a_y <-phenotyping$Growth_rate_mean[phenotyping$Genotype == "Anc"]
m_x <- phenotyping$Conjugation_rate_mean[phenotyping$Genotype == "Mut"]
m_y <-phenotyping$Growth_rate_mean[phenotyping$Genotype == "Mut"]

anc_inv <- readRDS("results/phenotyping/pB10_anc_inv.rds")
mut_inv <- readRDS("results/phenotyping/pB10_mut_inv.rds")

## Edits to plots ----
p1 <- mut_inv +
  labs(title = "Y invades X in LFC") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 10)) +
  guides(color = "none",
         shape = "none")

p2 <- LFC +
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.box.spacing = unit(0, "pt"),
        legend.spacing.y = unit(0, "pt")) +
  guides(color = guide_legend(order = 1),
         shape = guide_legend(order = 1))

p3 <- anc_inv +
  labs(title = "X invades Y in HFC") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 10)) +
  guides(color = "none",
         shape = "none")

p4 <- HFC +
  theme(legend.position = "none",
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(color = "black")) +
  guides(color = "none",
         shape = "none",
         linetype = "none")

## Final plot ----
final_plot <- p1 + free(p2) + p3 + p4 + 
  plot_layout(
    ncol = 2, 
    widths = c(1, 2)
  ) +
  plot_annotation(tag_levels = "A") & 
  theme(
    #plot.margin = margin(8, 8, 8, 8),
    plot.tag = element_text(size = 10)
  ) 

ggsave("figures/fig5_validation.pdf",
       final_plot, width = 5.75, height = 4.5, units = "in")

