# Generate Fig psweep for host differences

# Load packages ----
library(tidyverse)
library(patchwork)

# Read in figures ----
anc_HFC <- readRDS("results/parameter_sweeps/pHFC_S.pB10/pHFC_S.pB10_inv_change_plot2.rds")
mut_HFC <- readRDS("results/parameter_sweeps/pHFC_S.pB10-A/pHFC_S.pB10-A_inv_change_plot2.rds")
high_HFC <- readRDS("results/parameter_sweeps/pHFC_S.pB10-high/pHFC_S.pB10-high_inv_change_plot2.rds")
low_HFC <- readRDS("results/parameter_sweeps/pHFC_S.pB10-low/pHFC_S.pB10-low_inv_change_plot2.rds")


# Crop figures ----
## Determine range ----
ph <- read_csv("input_data/strain_phenotypes.csv")

# Select strains to plot
ph_hl <- data.frame("Strain" = c("S.pB10-high", "S.pB10-low"),
                    "Growth_rate" = c(0.39, 0.29),
                    "Conjugation_rate" = c(1.56e-12, 1.61e-14),
                    "Source" = c(NA, NA))
ph_plot <- ph %>%
  filter(Strain %in% c("S.pB10", "S.pB10-A")) %>%
  rbind(ph_hl) %>%
  mutate(log_Conj = log10(Conjugation_rate))

anc_gamma <- ph_plot$log_Conj[ph_plot$Strain == "S.pB10"]
mut_gamma <- ph_plot$log_Conj[ph_plot$Strain == "S.pB10-A"]
anc_psi <- ph_plot$Growth_rate[ph_plot$Strain == "S.pB10"]
mut_psi <- ph_plot$Growth_rate[ph_plot$Strain == "S.pB10-A"]

mid_gamma <- mean(c(anc_gamma, mut_gamma))
mid_psi <- mean(c(anc_psi, mut_psi))

max_gamma <- mid_gamma + 3
min_gamma <- mid_gamma - 3

max_psi <- mid_psi + 0.3
min_psi <- mid_psi - 0.3



## Crop figures, remove legends and labels from inside plots
fig_theme <- theme_bw() +
  theme(legend.position = "top")

pA <- mut_HFC +
  scale_x_continuous(limits = c(min_gamma, max_gamma)) +
  scale_y_continuous(limits = c(min_psi, max_psi)) +
  geom_point(data = ph_plot,
             mapping = aes(log_Conj, Growth_rate, shape = Strain),
             size = 3) +
  scale_shape_manual(values = c("S.pB10" = 16,
                                "S.pB10-A" = 17,
                                "S.pB10-high" = 0,
                                "S.pB10-low" = 5)) +
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        legend.title = element_blank())

pB <- high_HFC +
  scale_x_continuous(limits = c(min_gamma, max_gamma)) +
  scale_y_continuous(limits = c(min_psi, max_psi)) +
  geom_point(data = ph_plot,
             mapping = aes(log_Conj, Growth_rate, shape = Strain),
             size = 3) +
  scale_shape_manual(values = c("S.pB10" = 16,
                                "S.pB10-A" = 17,
                                "S.pB10-high" = 0,
                                "S.pB10-low" = 5)) +
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        legend.title = element_blank())

pC <- low_HFC +
  scale_x_continuous(limits = c(min_gamma, max_gamma)) +
  scale_y_continuous(limits = c(min_psi, max_psi)) +
  geom_point(data = ph_plot,
             mapping = aes(log_Conj, Growth_rate, shape = Strain),
             size = 3) +
  scale_shape_manual(values = c("S.pB10" = 16,
                                "S.pB10-A" = 17,
                                "S.pB10-high" = 0,
                                "S.pB10-low" = 5)) +
  theme(legend.title = element_blank())

pD <- anc_HFC +
  scale_x_continuous(limits = c(min_gamma, max_gamma)) +
  scale_y_continuous(limits = c(min_psi, max_psi)) +
  geom_point(data = ph_plot,
             mapping = aes(log_Conj, Growth_rate, shape = Strain),
             size = 3) +
  scale_shape_manual(values = c("S.pB10" = 16,
                                "S.pB10-A" = 17,
                                "S.pB10-high" = 0,
                                "S.pB10-low" = 5)) +
  theme(axis.title.y = element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_blank())


final_plot <- (pA | pB) /
  (pC | pD) +
  plot_layout(guides = "collect") &
  theme(legend.position = "top") 

ggsave("figures/fig_psweep_HFC_strain.pdf",
       final_plot, width = 10, height = 10, units = "in")



