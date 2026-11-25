# Process and plot HFC emulsion data

# Load packages ----
library(tidyverse)
library(broom)
library(ggpattern)
library(argparse)
library(jsonlite)

# Set arguments parser inputs ----
parser <- ArgumentParser()
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
invs_lw <- 0.25

## Points ----
plot_points <- jsonlite::fromJSON(args$points)
exp_point_size <- plot_points[["exp_point_size"]]
sh_Anc <- plot_points[["sh_Anc"]]
sh_Mut <- plot_points[["sh_Mut"]]

## Common aesthetics ----
anc_name <- "X"
mut_name <- "Y"

## Retrieve ggplot theme ----
source("src/ggplot_theme.R")

## Labels ----
axes_names <- c("Growth", "Conjugation", "Transconjugant\nselection")

# Read in data ----
## Specify output directory ----
output_dir <- "results/experimental_validation"
# Data processing/experimental parameters ----
tselect_t = 48 # emulsion grwoth time of 48 hr

## Relevel factors ----
host_levels <- c("rifR", "nalR")
ph_levels <- c("Anc", "Mut", "F")
phase_levels <- c("growth", "conj", "tselect")

# Read in data ----
plating_HFC <- read_csv("input_data/experimental_data/2026-10-19_HFC_emulsion_plating.csv")

extinction_out_av <- read_csv("results/experimental_validation/extinction_cell_counts/extinction_size_av.csv") %>%
  select(Host, Phenotype, Plate_type, Extinction_mean, Cells_colony_mean)

# Process HFC data ----
## Clean up data, apply extinctions, add metadata ----
plating_HFC2 <- plating_HFC %>%
  filter(!is.na(Count)) %>%
  mutate(Replicate = as.character(Replicate)) %>%
  mutate(Host = ifelse(str_detect(Plate_type, "Rif"), "rifR", "nalR")) %>% # Add host info
  left_join(extinction_out_av, by = c("Host", "Phenotype", "Plate_type")) %>% # Add extinction data
  mutate(Count_c = Count/(1-Extinction_mean)) %>% # Apply extinction
  mutate(Density_c = Count_c * (1000/Volume_plated) * (10^Plate_Dilution)) # Recalculate density

## Prep data for export ----
HFC_out <- plating_HFC2 %>%
  select(Experiment, Treatment, Phase, Time, Replicate, Phenotype, Host,
         Count, Count_c, Density_c, Count_type) %>%
  mutate(Phase = factor(Phase, levels = phase_levels),
         Host = factor(Host, levels = host_levels),
         Phenotype = factor(Phenotype, levels = ph_levels))

## Calculate frequencies ----
HFC_f <- HFC_out %>%
  filter(Phenotype %in% c("Anc", "Mut")) %>%
  filter(!(Time == 0 & Host == "nalR")) %>%
  filter(!(Time == 5 & Host == "rifR")) %>%
  select(-Count, -Count_c) %>%
  pivot_wider(names_from = Phenotype, values_from = Density_c) %>%
  mutate(Density_p = Anc + Mut,
         Anc_freq = Anc/Density_p,
         Mut_freq = Mut/Density_p) %>%
  pivot_longer(c(Anc_freq, Mut_freq), names_to = "Phenotype", values_to = "Frequency") %>%
  mutate(Phenotype = str_replace(Phenotype, "_freq", "")) %>%
  select(-Anc, -Mut)

## Average data ----
HFC_out_av <- HFC_out %>%
  group_by(Experiment, Treatment, Phase, Time, Phenotype, Host, Count_type) %>%
  summarize(Density_mean = mean(Density_c),
            Density_sd = sd(Density_c),
            n = n()) %>%
  ungroup() %>%
  mutate(Density_se = Density_sd/sqrt(n))

HFC_f_av <- HFC_f %>%
  group_by(Experiment, Treatment, Phase, Time, Phenotype, Host, Count_type) %>%
  summarize(Frequency_mean = mean(Frequency),
            Frequency_sd = sd(Frequency),
            n = n()) %>%
  ungroup() %>%
  mutate(Frequency_se = Frequency_sd/sqrt(n))

## Statistical analysis ----
HFC_stat <- HFC_f %>%
  filter(Time %in% c(0,5+tselect_t)) %>%
  filter(Phenotype == "Anc") %>%
  select(Time, Phenotype, Frequency, Replicate) %>%
  pivot_wider(names_from = Time, values_from = Frequency,
              names_prefix = "Freq_T") %>%
  mutate(Freq_T0a = asin(sqrt(Freq_T0)),
         !!paste0("Freq_T",5+tselect_t,"a") := asin(sqrt(!!sym(paste0("Freq_T",5+tselect_t))))) # arcsin transformation

HFC_tt <- t.test(HFC_stat[[paste0("Freq_T",5+tselect_t,"a")]], HFC_stat$Freq_T0a, 
                 alternative = "greater",
                 paired = TRUE,
                 var.equal = FALSE,
                 conf.level = 0.95)

## Save files ----
write_csv(HFC_out, paste(output_dir, "HFCe", "HFCe_plating_processed.csv", sep = "/"))
write_csv(HFC_out_av, paste(output_dir, "HFCe", "HFCe_plating_processed_av.csv", sep = "/"))
write_csv(HFC_f, paste(output_dir, "HFCe", "HFCe_frequency_processed.csv", sep = "/"))
write_csv(HFC_f_av, paste(output_dir, "HFCe", "HFCe_frequency_processed_av.csv", sep = "/"))
write_csv(tidy(HFC_tt), paste(output_dir, "HFCe", "HFCe_ttest.csv", sep = "/"))

# Plot data ----
phase_levels <- c("Growth", "Conjugation", "Transconjugant selection")
ph_levels <- c("Anc", "Mut", "F")

HFC_f_plot <- HFC_f_av %>%
  bind_rows(data.frame(Phase = rep("growth", 2),
                       Phenotype = c("Anc", "Mut"))) %>%
  mutate(Data = "experiment") %>%
  mutate(Phase_n = case_when(Phase == "growth" ~ 1,
                             Phase == "conj" & Time == 0 ~ 2,
                             Phase == "conj" & Time == 5 ~ 3,
                             Phase == "tselect" ~ 4)) %>%
  mutate(Phenotype = factor(Phenotype, levels = ph_levels))

HFC_f_plot2 <- HFC_f %>%
  mutate(Phase_n = case_when(Phase == "growth" ~ 1,
                             Phase == "conj" & Time == 0 ~ 2,
                             Phase == "conj" & Time == 5 ~ 3,
                             Phase == "tselect" ~ 4)) %>%
  mutate(Phenotype = factor(Phenotype, levels = ph_levels))

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
  geom_point(data = HFC_f_plot2,
             mapping = aes(Phase_n, Frequency, shape = Phenotype, color = Phenotype),
             size = exp_point_size + 0.5, alpha = 0.75) +
  geom_line(data = HFC_f_plot,
            mapping = aes(Phase_n, Frequency_mean, color = Phenotype, linetype = Data),
            linewidth = plot_lw) +
  geom_errorbar(data = HFC_f_plot,
                mapping = aes(x = Phase_n, ymax = Frequency_mean + Frequency_se, ymin = Frequency_mean - Frequency_se),
                width = 0.05, color = "black") +
  geom_point(data = filter(HFC_f_plot, Data == "experiment"),
             mapping = aes(Phase_n, Frequency_mean, shape = Phenotype),
             size = exp_point_size, color = "black") +
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
                                "Mut" = p_Mut),
                     labels = c("Anc" = anc_name,
                                "Mut" = mut_name)) +
  scale_shape_manual(values = c("Anc" = sh_Anc,
                                "Mut" = sh_Mut),
                     labels = c("Anc" = anc_name,
                                "Mut" = mut_name)) +
  scale_linetype_manual(values = c("simulation" = sim_l,
                                   "experiment" = exp_l)) +
  fig_aes +
  theme(axis.title.x = element_blank()) +
  guides(fill = "none")

ggsave(paste(output_dir, "HFCe", "HFCe_fig5d.pdf", sep = "/"), h1, 
       width = 4.5, height = 2.5, units = "in")
saveRDS(h1, paste(output_dir, "HFCe", "HFCe_fig5d.rds", sep = "/"))

h1_edit <- h1 +
  theme(legend.position = "top",
        legend.title = element_blank())

ggsave("figures/figS_HFC_emulsion.pdf", h1_edit, 
       width = 4.5, height = 3, units = "in")

