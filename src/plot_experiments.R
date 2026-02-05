# Plot experimental validation data

# Load packages ----
library(tidyverse)

# Set common aesthetics ----
p_Anc <- "#8394F6"
p_Mut <- "#8A407A"
p_F <- "gray40"

source("src/ggplot_theme.R")
p_conj <- "#D0D9FF"
p_tselect <- "#FBE9FF"
p_growth <- "#FFFBEF"

rifR_l <- "solid"
nalR_l <- "dashed"

# Read in data ----
## Specify output directory ----
output_dir <- "results/experimental_validation"

## Read files ----
HFC_out <- read_csv(paste(output_dir, "HFC", "HFC_plating_processed.csv", sep = "/"))
HFC_out_av <- read_csv(paste(output_dir, "HFC", "HFC_plating_processed_av.csv", sep = "/"))
HFC_f <- read_csv(paste(output_dir, "HFC", "HFC_frequency_processed.csv", sep = "/"))
HFC_f_av <- read_csv(paste(output_dir, "HFC", "HFC_frequency_processed_av.csv", sep = "/"))

LFC_out <- read_csv(paste(output_dir, "LFC", "LFC_plating_processed.csv", sep = "/"))
LFC_out_av <- read_csv(paste(output_dir, "LFC", "LFC_plating_processed_av.csv", sep = "/"))
LFC_f <- read_csv(paste(output_dir, "LFC", "LFC_frequency_processed.csv", sep = "/"))
LFC_f_av <- read_csv(paste(output_dir, "LFC", "LFC_frequency_processed_av.csv", sep = "/"))

# Plotting edits ----
## Replace 0s with 1 for easier log axis plotting ----
HFC_out$Density_c[HFC_out$Density_c == 0] <- 1
HFC_out_av$Density_mean[HFC_out_av$Density_mean == 0] <- 1
LFC_out$Density_c[LFC_out$Density_c == 0] <- 1
LFC_out_av$Density_mean[LFC_out_av$Density_mean == 0] <- 1

## Relevel factors ----
phase_levels <- c("Growth", "Conjugation", "Transconjugant selection")

## Generate datasets for background plots ----
HFC_phases <- data.frame("Phase" = c("Conjugation","Transconjugant selection"),
                         "T_start" = c(0,5),
                         "T_end" = c(5,77)) %>%
  mutate(Phase = factor(Phase, levels = phase_levels))

LFC_phases <- data.frame("Phase" = c("Growth","Conjugation","Transconjugant selection"),
                         "T_start" = c(0,168,173),
                         "T_end" = c(168,173,245)) %>%
  mutate(Phase = factor(Phase, levels = phase_levels))


# Plots ----
## HFC ----

### Density plot ----
h1 <- ggplot() +
  geom_rect(data = HFC_phases,
            mapping = aes(xmin = T_start, xmax = T_end, ymin = 0, ymax = Inf, fill = Phase)) +
  scale_fill_manual(values = c("Conjugation" = p_conj,
                               "Transconjugant selection" = p_tselect)) +
  geom_line(data = HFC_out_av,
            mapping = aes(Time, Density_mean, color = Phenotype, linetype = Host),
            linewidth = 2) +
  geom_errorbar(data = HFC_out_av,
                mapping = aes(x = Time, ymax = Density_mean + Density_se, ymin = Density_mean - Density_se,
                              color = Phenotype, linetype = Host),
                width = 1) +
  geom_point(data = HFC_out_av,
             mapping = aes(Time, Density_mean, color = Phenotype, shape = Count_type),
             size = 4) +
  scale_shape_manual(values = c("count" = 16,
                                "estimate" = 15,
                                "below_limit" = 17,
                                "dilution_est" = 18)) +
  scale_y_continuous(trans = "log10", name = "Density (CFU/mL)",
                     breaks = c(1e0,1e2,1e4,1e6,1e8),
                     labels = sapply(c(0,2,4,6,8),function(i){parse(text = sprintf("10^%d",i))})) +
  scale_x_continuous(expand = c(0.001, 0.001)) + 
  scale_color_manual(values = c("Anc" = p_Anc,
                                "Mut" = p_Mut,
                                "F" = p_F)) +
  scale_linetype_manual(values = c("rifR" = rifR_l,
                                   "nalR" = nalR_l)) +
  fig_aes

### Frequency plot ----
h2 <- ggplot() +
  geom_rect(data = HFC_phases,
            mapping = aes(xmin = T_start, xmax = T_end, ymin = 0, ymax = Inf, fill = Phase)) +
  scale_fill_manual(values = c("Conjugation" = p_conj,
                               "Transconjugant selection" = p_tselect)) +
  geom_line(data = HFC_f_av,
            mapping = aes(Time, Frequency_mean, color = Phenotype),
            linewidth = 2) +
  geom_errorbar(data = HFC_f_av,
                mapping = aes(x = Time, ymax = Frequency_mean + Frequency_se, ymin = Frequency_mean - Frequency_se,
                              color = Phenotype),
                width = 1) +
  geom_point(data = HFC_f_av,
             mapping = aes(Time, Frequency_mean, color = Phenotype, shape = Count_type),
             size = 4) +
  scale_shape_manual(values = c("count" = 16,
                                "estimate" = 15,
                                "below_limit" = 17,
                                "dilution_est" = 18)) +
  scale_y_continuous(limits = c(0.00001,1),
                     breaks=c(0, 0.0001, 0.01, 1),
                     labels=c(0,sapply(c(-4,-2),function(i){parse(text = sprintf("10^%d",i))}),1),
                     expand = c(0.01, 0.01),
                     trans = "log10", name = "Frequency") +
  scale_x_continuous(expand = c(0.001, 0.001)) +
  scale_color_manual(values = c("Anc" = p_Anc,
                                "Mut" = p_Mut)) +
  fig_aes

## LFC ----

### Density plot ----
l1 <- ggplot() +
  geom_rect(data = LFC_phases,
            mapping = aes(xmin = T_start, xmax = T_end, ymin = 0, ymax = Inf, fill = Phase)) +
  scale_fill_manual(values = c("Conjugation" = p_conj,
                               "Transconjugant selection" = p_tselect,
                               "Growth" = p_growth)) +
  geom_line(data = LFC_out_av,
            mapping = aes(Time, Density_mean, color = Phenotype, linetype = Host),
            linewidth = 2) +
  geom_errorbar(data = LFC_out_av,
                mapping = aes(x = Time, ymax = Density_mean + Density_se, ymin = Density_mean - Density_se,
                              color = Phenotype, linetype = Host),
                width = 1) +
  geom_point(data = LFC_out_av,
             mapping = aes(Time, Density_mean, color = Phenotype, shape = Count_type),
             size = 4) +
  scale_shape_manual(values = c("count" = 16,
                                "estimate" = 15,
                                "below_limit" = 17,
                                "dilution_est" = 18)) +
  scale_y_continuous(trans = "log10", name = "Density (CFU/mL)",
                     breaks = c(1e0,1e2,1e4,1e6,1e8),
                     labels = sapply(c(0,2,4,6,8),function(i){parse(text = sprintf("10^%d",i))})) +
  scale_x_continuous(expand = c(0.001, 0.001)) +
  scale_color_manual(values = c("Anc" = p_Anc,
                                "Mut" = p_Mut,
                                "F" = p_F)) +
  scale_linetype_manual(values = c("rifR" = rifR_l,
                                   "nalR" = nalR_l)) +
  fig_aes

### Frequency plot ----
l2 <- ggplot() +
  geom_rect(data = LFC_phases,
            mapping = aes(xmin = T_start, xmax = T_end, ymin = 0, ymax = Inf, fill = Phase)) +
  scale_fill_manual(values = c("Conjugation" = p_conj,
                               "Transconjugant selection" = p_tselect,
                               "Growth" = p_growth)) +
  geom_line(data = LFC_f_av,
            mapping = aes(Time, Frequency_mean, color = Phenotype),
            linewidth = 2) +
  geom_errorbar(data = LFC_f_av,
                mapping = aes(x = Time, ymax = Frequency_mean + Frequency_se, ymin = Frequency_mean - Frequency_se,
                              color = Phenotype),
                width = 1) +
  geom_point(data = LFC_f_av,
             mapping = aes(Time, Frequency_mean, color = Phenotype, shape = Count_type),
             size = 4) +
  scale_shape_manual(values = c("count" = 16,
                                "estimate" = 15,
                                "below_limit" = 17,
                                "dilution_est" = 18)) +
  scale_y_continuous(limits = c(0.00001,1),
                     breaks=c(0, 0.0001, 0.01, 1),
                     labels=c(0,sapply(c(-4,-2),function(i){parse(text = sprintf("10^%d",i))}),1),
                     expand = c(0.01, 0.01),
                     trans = "log10", name = "Frequency") +
  scale_x_continuous(expand = c(0.001, 0.001)) +
  scale_color_manual(values = c("Anc" = p_Anc,
                                "Mut" = p_Mut)) +
  fig_aes

# Save plots ----
ggsave(paste(output_dir, "HFC", "HFC_density_plot.pdf", sep = "/"), h1, 
       width = 12, height = 6, units = "in")
ggsave(paste(output_dir, "HFC", "HFC_frequency_plot.pdf", sep = "/"), h2, 
       width = 12, height = 6, units = "in")

ggsave(paste(output_dir, "LFC", "LFC_density_plot.pdf", sep = "/"), l1, 
       width = 16, height = 6, units = "in")
ggsave(paste(output_dir, "LFC", "LFC_frequency_plot.pdf", sep = "/"), l2, 
       width = 16, height = 6, units = "in")

