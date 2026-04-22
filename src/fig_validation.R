# Generate experimental validation figure

# Load packages ----
library(tidyverse)
library(patchwork)

# Read in figures ----
HFC_sim <- readRDS("results/case_study_sims/HFC_S.pB10/HFC_S.pB10_frequency_plot.rds")
LFC_sim <- readRDS("results/case_study_sims/LFC_S.pB10-A/LFC_S.pB10-A_frequency_plot.rds")
HFC_v <- readRDS("results/experimental_validation/HFC/HFC_frequency_plot.rds")
LFC_v <- readRDS("results/experimental_validation/LFC/LFC_frequency_plot.rds")

# Edits to plots ----
pA <- LFC_sim +
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "top",
        legend.title = element_blank())

pB <- HFC_sim +
  guides(shape = "none",
         fill = "none",
         color = "none") +
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none")

pC <- LFC_v +
  labs(x = "Time (hr)",
       y = "Frequency",
       color = "Genotype",
       fill = "Phase") +
  guides(shape = "none",
         fill = "none",
         color = "none") +
  theme(legend.position = "none")

pD <- HFC_v +
  labs(x = "Time (hr)",
       y = "Frequency",
       color = "Genotype",
       fill = "Phase") +
  guides(shape = "none",
         fill = "none",
         color = "none") +
  theme(axis.title.y = element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none")

# Final plot
top_row    <- pA | pB + plot_layout(widths = c(2, 1))
bottom_row <- pC | pD + plot_layout(widths = c(2, 1))

final_plot <- (pA / pC | 
               pB / pD) +
  plot_layout(ncol = 2,
              widths = c(2,1),
              guides = "collect",
              axes = "collect",
              axis_titles = "collect") &
  theme(legend.position = "top") 

ggsave("figures/fig5_validation.pdf",
       final_plot, width = 14, height = 8, units = "in")



# Second version with invasion comparison plots
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

a_x <- phenotyping$Conjugation_rate_mean[phenotyping$Genotype == "Anc"]
a_y <-phenotyping$Growth_rate_mean[phenotyping$Genotype == "Anc"]
m_x <- phenotyping$Conjugation_rate_mean[phenotyping$Genotype == "Mut"]
m_y <-phenotyping$Growth_rate_mean[phenotyping$Genotype == "Mut"]

anc_inv <- readRDS("results/phenotyping/pB10_anc_inv.rds")
mut_inv <- readRDS("results/phenotyping/pB10_mut_inv.rds")

# Edits to plots ----
p1 <- anc_inv +
  theme(legend.position = "none") +
  guides(color = "none")

p2 <- HFC_sim +
  guides(shape = "none",
         fill = "none",
         color = "none") +
  theme(axis.title.x = element_blank(),
        legend.position = "none")

p3 <- HFC_v +
  labs(x = "Time (hr)",
       y = "Frequency",
       color = "Genotype",
       fill = "Phase") +
  guides(shape = "none",
         fill = "none",
         color = "none") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none")

p4 <- mut_inv +
  theme(legend.position = "none") +
  guides(color = "none")

p5 <- LFC_sim +
  theme(legend.position = "bottom",
        legend.title = element_blank())

p6 <- LFC_v +
  labs(x = "Time (hr)",
       y = "Frequency",
       color = "Genotype",
       fill = "Phase") +
  guides(shape = "none",
         fill = "none",
         color = "none") +
  theme(axis.title.y = element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none")


# Final plot
final_plot <- p1 + p2 + p3 + p4 + p5 + p6 + 
  plot_layout(
    ncol = 3, 
    widths = c(1, 2, 2),
    guides = "collect"
  ) & 
  theme(
    legend.position = "top",
    legend.box = "horizontal",
    plot.margin = margin(10, 10, 10, 10)
  ) 

ggsave("figures/fig5_validation2.pdf",
       final_plot, width = 19.5, height = 9, units = "in")

