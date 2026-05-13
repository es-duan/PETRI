# Plot experimental validation data

# Load packages ----
library(tidyverse)
library(argparse)
library(jsonlite)

# Set arguments parser inputs ----
parser <- ArgumentParser()
parser$add_argument("-t","--treatment", type = "character", help = "Specify Treatment ID")
parser$add_argument("-c","--colors", help = "JSON string of plot colors")
parser$add_argument("-l","--lines", help = "JSON string of plot lines")

# Parse arguments
args <- parser$parse_args()

# Load global variables ----
## Colors ----
plot_colors <- fromJSON(args$colors)
p_Anc <- plot_colors[["p_Anc"]]
p_Mut <- plot_colors[["p_Mut"]]
p_F <- plot_colors[["p_F"]]
p_growth <- plot_colors[["p_growth"]]
p_conj <- plot_colors[["p_conj"]]
p_tselect <- plot_colors[["p_tselect"]]
p_imm <- plot_colors[["p_imm"]]

## Lines ----
plot_lines <- fromJSON(args$lines)
rifR_l <- plot_lines[["rifR_l"]]
nalR_l <- plot_lines[["nalR_l"]]

## Retrieve ggplot theme ----
source("src/ggplot_theme.R")

# Read in data ----
## Specify output directory ----
output_dir <- "results/experimental_validation"

## Read files ----
HFC_out_av <- read_csv(paste(output_dir, "HFC", "HFC_plating_processed_av.csv", sep = "/"))
HFC_f_av <- read_csv(paste(output_dir, "HFC", "HFC_frequency_processed_av.csv", sep = "/"))

LFC_out_av <- read_csv(paste(output_dir, "LFC", "LFC_plating_processed_av.csv", sep = "/"))
LFC_f_av <- read_csv(paste(output_dir, "LFC", "LFC_frequency_processed_av.csv", sep = "/"))

# Plotting edits ----
## Replace 0s with 1 for easier log axis plotting ----
HFC_out_av$Density_mean[HFC_out_av$Density_mean == 0] <- 1
LFC_out_av$Density_mean[LFC_out_av$Density_mean == 0] <- 1

## Relevel factors ----
phase_levels <- c("Growth", "Conjugation", "Transconjugant selection")
ph_levels <- c("Anc", "Mut", "F")

## Generate datasets for background plots ----
tselect_t <- 72
HFC_phases <- data.frame("Phase" = c("Conjugation","Transconjugant selection"),
                         "T_start" = c(0,5),
                         "T_end" = c(5,5 + tselect_t)) %>%
  mutate(Phase = factor(Phase, levels = phase_levels))

LFC_phases <- data.frame("Phase" = c("Growth","Conjugation","Transconjugant selection"),
                         "T_start" = c(0,168,173),
                         "T_end" = c(168,173,173 + tselect_t)) %>%
  mutate(Phase = factor(Phase, levels = phase_levels))

## Divide datasets to differentiate strains to be selected against ----
HFC_out_av2 <- HFC_out_av %>%
  mutate(selection = case_when(Host == "nalR" & Phenotype %in% c("Anc", "Mut") ~ "no",
                               TRUE ~ "yes")) %>%
  mutate(Phenotype = factor(Phenotype, levels = ph_levels))
HFC_B <- HFC_out_av2 %>% filter(selection == "no")
HFC_L <- HFC_out_av2 %>% filter(selection == "yes")
HFC_P <- HFC_L %>% filter(Time == 5)
HFC_P2 <- HFC_L %>% filter(Time != 5) 

LFC_out_av2 <- LFC_out_av %>%
  mutate(selection = case_when(Host == "nalR" & Phenotype %in% c("Anc", "Mut") ~ "no",
                               Phase == "growth" ~ "no",
                               TRUE ~ "yes")) %>%
  mutate(Phenotype = factor(Phenotype, levels = ph_levels))
LFC_L <- LFC_out_av2 %>% filter(selection == "yes")
LFC_B <- LFC_out_av2 %>% 
  filter(selection == "yes" & Time == 168) %>%
  mutate(Time = ifelse(selection == "yes" & Time == 168, 168.001, Time)) %>%
  rbind(filter(LFC_out_av2, selection == "no"))
LFC_P <- LFC_L %>% filter(Time == 173)
LFC_P2 <- LFC_L %>% filter(Time != 173)


# Plots ----
## HFC ----

### Density plot ----
h1 <- ggplot() +
  geom_rect(data = HFC_phases,
            mapping = aes(xmin = T_start, xmax = T_end, ymin = 0, ymax = Inf, fill = Phase)) +
  scale_fill_manual(values = c("Conjugation" = p_conj,
                               "Transconjugant selection" = p_tselect)) +
  geom_line(data = HFC_L,
            mapping = aes(Time, Density_mean, color = Phenotype, linetype = Host),
            linewidth = 2, alpha = 0.5) +
  geom_errorbar(data = HFC_L,
                mapping = aes(x = Time, ymax = Density_mean + Density_se, ymin = Density_mean - Density_se,
                              color = Phenotype, linetype = Host),
                width = 1, alpha = 0.5) +
  geom_point(data = HFC_P2,
             mapping = aes(Time, Density_mean, color = Phenotype),
             size = 4, alpha = 0.75) +
  geom_point(data = HFC_P,
             mapping = aes(Time, Density_mean, color = Phenotype),
             size = 4, shape = 21, fill = "white", alpha = 0.75) +
  geom_line(data = HFC_B,
            mapping = aes(Time, Density_mean, color = Phenotype, linetype = Host),
            linewidth = 2) +
  geom_errorbar(data = HFC_B,
                mapping = aes(x = Time, ymax = Density_mean + Density_se, ymin = Density_mean - Density_se,
                              color = Phenotype, linetype = Host),
                width = 0.3) +
  geom_point(data = HFC_B,
             mapping = aes(Time, Density_mean, color = Phenotype, shape = Count_type),
             size = 4) +
  scale_shape_manual(values = c("count" = 16,
                                "estimate" = 15,
                                "below_limit" = 17,
                                "dilution_est" = 18)) +
  scale_y_continuous(trans = "log10", name = "Density (CFU/mL)",
                     breaks = c(1e0,1e2,1e4,1e6,1e8),
                     labels = sapply(c(0,2,4,6,8),function(i){parse(text = sprintf("10^%d",i))}),
                     limits = c(0.9, 1.5e9),
                     expand = c(0.01, 0.01)) +
  scale_x_continuous(expand = c(0.001, 0.001)) + 
  scale_color_manual(values = c("Anc" = p_Anc,
                                "Mut" = p_Mut,
                                "F" = p_F)) +
  scale_linetype_manual(values = c("rifR" = rifR_l,
                                   "nalR" = nalR_l)) +
  guides(color = guide_legend(order=1, override.aes = list(shape = NA))) +
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
                width = 0.3) +
  geom_point(data = HFC_f_av,
             mapping = aes(Time, Frequency_mean, color = Phenotype),
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
  geom_line(data = LFC_L,
            mapping = aes(Time, Density_mean, color = Phenotype, linetype = Host),
            linewidth = 2, alpha = 0.5) +
  geom_errorbar(data = LFC_L,
                mapping = aes(x = Time, ymax = Density_mean + Density_se, ymin = Density_mean - Density_se,
                              color = Phenotype, linetype = Host),
                width = 1, alpha = 0.5) +
  geom_point(data = LFC_P,
             mapping = aes(Time, Density_mean, color = Phenotype, shape = Count_type),
             size = 4, alpha = 0.75) +
  geom_point(data = LFC_P2,
             mapping = aes(Time, Density_mean, color = Phenotype, shape = Count_type),
             size = 4, alpha = 0.75) +
  geom_line(data = LFC_B,
            mapping = aes(Time, Density_mean, color = Phenotype, linetype = Host),
            linewidth = 2) +
  geom_errorbar(data = LFC_B,
                mapping = aes(x = Time, ymax = Density_mean + Density_se, ymin = Density_mean - Density_se,
                              color = Phenotype, linetype = Host),
                width = 1) +
  geom_point(data = LFC_B,
             mapping = aes(Time, Density_mean, color = Phenotype, shape = Count_type),
             size = 4) +
  scale_shape_manual(values = c("count" = 16,
                                "estimate" = 15,
                                "below_limit" = 17,
                                "dilution_est" = 18)) +
  scale_y_continuous(trans = "log10", name = "Density (CFU/mL)",
                     breaks = c(1e0,1e2,1e4,1e6,1e8),
                     labels = sapply(c(0,2,4,6,8),function(i){parse(text = sprintf("10^%d",i))}),
                     limits = c(0.9, 1.5e9),
                     expand = c(0.01, 0.01)) +
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
             mapping = aes(Time, Frequency_mean, color = Phenotype),
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
saveRDS(h1, paste(output_dir, "HFC", "HFC_density_plot.rds", sep = "/"))
ggsave(paste(output_dir, "HFC", "HFC_frequency_plot.pdf", sep = "/"), h2, 
       width = 12, height = 6, units = "in")
saveRDS(h2, paste(output_dir, "HFC", "HFC_frequency_plot.rds", sep = "/"))


ggsave(paste(output_dir, "LFC", "LFC_density_plot.pdf", sep = "/"), l1, 
       width = 16, height = 6, units = "in")
saveRDS(l1, paste(output_dir, "LFC", "LFC_density_plot.rds", sep = "/"))
ggsave(paste(output_dir, "LFC", "LFC_frequency_plot.pdf", sep = "/"), l2, 
       width = 16, height = 6, units = "in")
saveRDS(l2, paste(output_dir, "LFC", "LFC_frequency_plot.rds", sep = "/"))

