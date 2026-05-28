# Plot experimental validation frequencies by phase

# Load packages ----
library(tidyverse)
library(ggpattern)
library(argparse)
library(jsonlite)

# Set arguments parser inputs ----
parser <- ArgumentParser()
parser$add_argument("-t","--treatment", type = "character", help = "Specify Treatment ID")
parser$add_argument("-c","--colors", help = "JSON string of plot colors")
parser$add_argument("-l","--lines", help = "JSON string of plot lines")
parser$add_argument("-o","--points", help = "JSON string of point aesthetics")

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
exp_l <- plot_lines[["exp_l"]]
sim_l <- plot_lines[["sim_l"]]
plot_lw <- plot_lines[["plot_lw"]]

## Points ----
plot_points <- jsonlite::fromJSON(args$points)
exp_point_size <- plot_points[["exp_point_size"]]
sh_Anc <- plot_points[["sh_Anc"]]
sh_Mut <- plot_points[["sh_Mut"]]

## Retrieve ggplot theme ----
source("src/ggplot_theme.R")

## Labels ----
axes_names <- c("Growth", "Conjugation", "Transconjugant\nselection")

# Read in data ----
## Read files ----
HFC_f_av <- read_csv("results/experimental_validation/HFC/HFC_frequency_processed_av.csv")
LFC_f_av <- read_csv("results/experimental_validation/LFC/LFC_frequency_processed_av.csv")

HFC_tt <- read_csv("results/experimental_validation/HFC/HFC_ttest.csv")
LFC_tt <- read_csv("results/experimental_validation/LFC/LFC_ttest.csv")

HFC_f_sim <- read_csv("results/case_study_sims/HFC_S.pB10/HFC_S.pB10_frequency_plot_df.csv")
LFC_f_sim <- read_csv("results/case_study_sims/LFC_S.pB10-A/LFC_S.pB10-A_frequency_plot_df.csv")

## Relevel factors ----
phase_levels <- c("Growth", "Conjugation", "Transconjugant selection")
ph_levels <- c("Anc", "Mut", "F")

# Select the first and last point for each phase ----
## HFC: add growth phase ----
# Merge with sim data
HFC_f_sim2 <- HFC_f_sim %>%
  filter(!(Phase == "tselect" & Time == 5)) %>%
  rename(Frequency_mean = Frequency) %>%
  mutate(Phenotype = ifelse(Genotype == "Ancestor", "Anc", "Mut")) %>%
  select(-Cycle, -Genotype) %>%
  bind_rows(data.frame(Phase = rep("growth", 2),
                       Phenotype = c("Anc", "Mut"))) %>%
  mutate(Data = "simulation")

HFC_f_plot <- HFC_f_av %>%
  bind_rows(data.frame(Phase = rep("growth", 2),
                       Phenotype = c("Anc", "Mut"))) %>%
  mutate(Data = "experiment") %>%
  bind_rows(HFC_f_sim2) %>%
  mutate(Phase_n = case_when(Phase == "growth" ~ 1,
                             Phase == "conj" & Time == 0 ~ 2,
                             Phase == "conj" & Time == 5 ~ 3,
                             Phase == "tselect" ~ 4)) %>%
  mutate(Phenotype = factor(Phenotype, levels = ph_levels))

hfc_s <- HFC_f_plot$Frequency_mean[HFC_f_plot$Time == 0 & HFC_f_plot$Phenotype == "Anc"][1]
hfc_e <- HFC_f_plot$Frequency_mean[HFC_f_plot$Time == 29 & HFC_f_plot$Phenotype == "Anc"][1]

## LFC: select first and last point of growth phase ----
# Merge with sim data
LFC_f_sim2 <- LFC_f_sim %>%
  filter(Phase == "growth" & Time %in% c(0, 168) |
           Phase == "conj" & Time == 173 |
           Phase == "tselect" & Time == 245) %>%
  rename(Frequency_mean = Frequency) %>%
  mutate(Phenotype = ifelse(Genotype == "Ancestor", "Anc", "Mut")) %>%
  select(-Cycle, -Genotype) %>%
  mutate(Data = "simulation")

LFC_f_plot <- LFC_f_av %>%
  filter(Phase != "growth" | Phase == "growth" & Time %in% c(0, 168)) %>%
  mutate(Data = "experiment") %>%
  bind_rows(LFC_f_sim2) %>%
  mutate(Phase_n = case_when(Phase == "growth" & Time == 0 ~ 1,
                             Phase == "growth" & Time == 168 ~ 2,
                             Phase == "conj" ~ 3,
                             Phase == "tselect" ~ 4)) %>%
  mutate(Phenotype = factor(Phenotype, levels = ph_levels))

lfc_s <- LFC_f_plot$Frequency_mean[LFC_f_plot$Time == 0 & LFC_f_plot$Phenotype == "Mut"][1]
lfc_e <- LFC_f_plot$Frequency_mean[LFC_f_plot$Time == 197 & LFC_f_plot$Phenotype == "Mut"][1]

# Datasets for background plots ----
phases <- data.frame("Phase" = c("Growth","Conjugation","Transconjugant selection"),
                     "T_start" = c(-Inf,2,3),
                     "T_end" = c(2,3,Inf)) %>%
  mutate(Phase = factor(Phase, levels = phase_levels))

# Plot frequencies by phase ----
## HFC ----
h1 <- ggplot() +
  geom_rect(data = phases,
            mapping = aes(xmin = T_start, xmax = T_end, ymin = 0, ymax = Inf, fill = Phase)) +
  geom_rect_pattern(aes(xmin = -Inf, xmax = 2, ymin = 0, ymax = Inf),
                    pattern = "stripe", pattern_color = "gray60", fill = NA, pattern_density = 0.09,
                    pattern_spacing = 0.03) +
  scale_fill_manual(values = c("Conjugation" = p_conj,
                               "Transconjugant selection" = p_tselect,
                               "Growth" = p_growth)) +
  # annotate("segment",
  #          x = 2, xend = 4, y = hfc_s, yend = hfc_s,
  #          color = "white", linewidth = 1) +
  # annotate("segment",
  #          x = 4, xend = 4, y = hfc_s, yend = hfc_e,
  #          color = "white", linewidth = 1) +
  # annotate("text", label = paste0("p = ", round(HFC_tt$p.value, 2)),
  #          x = 3, hfc_s - 0.001,
  #          color = "white") +
  geom_line(data = HFC_f_plot,
            mapping = aes(Phase_n, Frequency_mean, color = Phenotype, linetype = Data),
            linewidth = plot_lw) +
  geom_errorbar(data = HFC_f_plot,
                mapping = aes(x = Phase_n, ymax = Frequency_mean + Frequency_se, ymin = Frequency_mean - Frequency_se,
                              color = Phenotype),
                width = 0.05) +
  geom_point(data = filter(HFC_f_plot, Data == "experiment"),
             mapping = aes(Phase_n, Frequency_mean, color = Phenotype, shape = Phenotype),
             size = exp_point_size) +
  scale_y_continuous(limits = c(0.00001,1),
                     breaks=c(0, 0.0001, 0.01, 1),
                     labels=c(0,sapply(c(-4,-2),function(i){parse(text = sprintf("10^%d",i))}),1),
                     expand = c(0.05, 0.05),
                     trans = "log10", name = "Frequency") +
  scale_x_continuous(expand = c(0.02, 0.02),
                     breaks = c(1.5, 2.5, 3.5),
                     labels = axes_names,
                     limits = c(1, 4)) +
  scale_color_manual(values = c("Anc" = p_Anc,
                                "Mut" = p_Mut)) +
  scale_linetype_manual(values = c("simulation" = sim_l,
                                   "experiment" = exp_l)) +
  fig_aes +
  theme(axis.title.x = element_blank()) +
  guides(fill = "none")

ggsave("figures/panels/fig5d_HFC.pdf", h1, 
       width = 4.5, height = 2.5, units = "in")
saveRDS(h1, "figures/panels/fig5d_HFC.rds")

## LFC ----
l1 <- ggplot() +
  geom_rect(data = phases,
            mapping = aes(xmin = T_start, xmax = T_end, ymin = 0, ymax = Inf, fill = Phase)) +
  scale_fill_manual(values = c("Conjugation" = p_conj,
                               "Transconjugant selection" = p_tselect,
                               "Growth" = p_growth)) +
  # annotate("segment",
  #          x = 1, xend = 4, y = lfc_s, yend = lfc_s,
  #          color = "white", linewidth = 1) +
  # annotate("segment",
  #          x = 4, xend = 4, y = lfc_s, yend = lfc_e,
  #          color = "white", linewidth = 1) +
  # annotate("text", label = "p < 0.01",
  #          x = 3, hfc_s - 0.001,
  #          color = "white") +
  geom_line(data = LFC_f_plot,
            mapping = aes(Phase_n, Frequency_mean, color = Phenotype, linetype = Data),
            linewidth = plot_lw) +
  geom_errorbar(data = LFC_f_plot,
                mapping = aes(x = Phase_n, ymax = Frequency_mean + Frequency_se, ymin = Frequency_mean - Frequency_se,
                              color = Phenotype),
                width = 0.05) +
  geom_point(data = filter(LFC_f_plot, Data == "experiment"),
             mapping = aes(Phase_n, Frequency_mean, color = Phenotype, shape = Phenotype),
             size = exp_point_size) +
  scale_y_continuous(limits = c(0.00001,1),
                     breaks=c(0, 0.0001, 0.01, 1),
                     labels=c(0,sapply(c(-4,-2),function(i){parse(text = sprintf("10^%d",i))}),1),
                     expand = c(0.05, 0.05),
                     trans = "log10", name = "Frequency") +
  scale_x_continuous(expand = c(0.02, 0.02),
                     breaks = c(1.5, 2.5, 3.5),
                     labels = axes_names,
                     limits = c(1,4)) +
  scale_color_manual(values = c("Anc" = p_Anc,
                                "Mut" = p_Mut)) +
  scale_linetype_manual(values = c("simulation" = sim_l,
                                   "experiment" = exp_l)) +
  fig_aes +
  theme(axis.title.x = element_blank()) +
  guides(fill = "none")

ggsave("figures/panels/fig5b_LFC.pdf", l1, 
       width = 4.5, height = 2.5, units = "in")
saveRDS(l1, "figures/panels/fig5b_LFC.rds")

