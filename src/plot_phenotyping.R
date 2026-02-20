# Plot phenotyping data

# Load packages ----
library(tidyverse)


# Load in data ----
LDM_phenotype <- read_csv("results/phenotyping/LDM_phenotyping.csv")
LDM_phenotype_av <- read_csv("results/phenotyping/LDM_phenotyping_av.csv")


# Set common aesthetics ----
p_Anc <- "#8394F6"
p_Mut <- "#8A407A"

source("src/ggplot_theme.R")

# Plot growth and conjugation rates ----
p1 <- ggplot() +
  geom_errorbar(data = LDM_phenotype_av,
                mapping= aes(x = Conjugation_rate_mean,
                             ymin = Growth_rate_mean - Growth_rate_se,
                             ymax = Growth_rate_mean + Growth_rate_se,
                             color = Genotype),
                width = 0) +
  geom_errorbarh(data = LDM_phenotype_av,
                 mapping = aes(y = Growth_rate_mean,
                               xmin = Conjugation_rate_mean - Conjugation_rate_se,
                               xmax = Conjugation_rate_mean + Conjugation_rate_se,
                               color = Genotype),
                 height = 0) +
  geom_point(data = LDM_phenotype_av,
             mapping = aes(Conjugation_rate_mean, Growth_rate_mean, color = Genotype),
             size = 3) +
  scale_x_continuous(trans = "log10", 
                     name = "Conjugation rate",
                     limits = c(1e-15, 1e-11),
                     breaks = c(1e-15,1e-13,1e-11),
                     labels = sapply(c(-15,-13,-11),function(i){parse(text = sprintf("10^%d",i))})) +
  scale_y_continuous(name = "Growth Rate",
                     limits = c(0.1, 0.5)) + 
  scale_color_manual(values = c("Anc" = p_Anc,
                                "Mut" = p_Mut)) +
  fig_aes

ggsave("results/phenotyping/pB10_phenotyping.pdf",
       p1, width = 8, height = 6, units = "in")
