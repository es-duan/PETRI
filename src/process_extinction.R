# Process and plot Extinction and Colony Size data

# Load packages ----
library(tidyverse)
library(argparse)
library(jsonlite)

# Set arguments parser inputs ----
parser <- ArgumentParser()
parser$add_argument("-t","--treatment", type = "character", help = "Specify Treatment ID")
parser$add_argument("-c","--colors", help = "JSON string of plot colors")
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

## Points ----
plot_points <- jsonlite::fromJSON(args$points)
exp_point_size <- plot_points[["exp_point_size"]]
sh_Anc <- plot_points[["sh_Anc"]]
sh_Mut <- plot_points[["sh_Mut"]]
sh_F <- plot_points[["sh_F"]]

## Retrieve ggplot theme ----
source("src/ggplot_theme.R")

# Read in data ----
plating_T <- read_csv("input_data/experimental_data/2026-01-27_T_extinction_plating.csv")
plating_DR <- read.csv("input_data/experimental_data/2026-01-27_DR_extinction_plating.csv")

# Specify output directory
output_dir <- "results/experimental_validation/extinction_cell_counts"

# T: Calculate colony size and extinctions ----
size <- plating_T %>%
  filter(Plate_type == "Tet50/Nal50") %>%
  mutate(Sizes = paste(Colony_size, Plated_size, sep = "-")) %>%
  mutate(Cells_colony = Count * (1000/Volume_plated) * 10^(Plate_Dilution) * Volume_suspended)

extinction <- plating_T %>%
  mutate(Sizes = paste(Colony_size, Plated_size, sep = "-")) %>%
  select(Experiment, Phenotype, Treatment, Sizes, Replicate, Plate_type, Count) %>%
  pivot_wider(names_from = Plate_type, values_from = Count) %>%
  mutate(Extinction = 1 - (`Tet50/Nal50`/LB)) %>%
  # Replace negative values with 0
  mutate(Extinction = ifelse(Extinction < 0, 0, Extinction))

processed <- size %>%
  left_join(extinction, by = c("Experiment", "Phenotype", "Treatment", "Sizes", "Replicate"))

# Reformat data for plotting ----
extinction_T <- processed %>%
  select(Experiment, Treatment, Phenotype, Sizes,
         Replicate, Cells_colony, Extinction) %>%
  mutate(Host = "nalR") %>%
  mutate(Plate_type = "Tet50/Nal50")

extinction_T_av <- extinction_T %>%
  group_by(Host, Phenotype, Plate_type) %>%
  summarise(Cells_colony_mean = mean(Cells_colony),
            Cells_colony_sd = sd(Cells_colony),
            Extinction_mean = mean(Extinction),
            Extinction_sd = sd(Extinction),
            n = n()) %>%
  ungroup() %>%
  mutate(Cells_colony_se = Cells_colony_sd/sqrt(n),
         Extinction_se = Extinction_sd/sqrt(n))

## Plots with final data ----
### Size plot ----
s3 <- ggplot() +
  # geom_point(data = extinction_T, 
  #            mapping = aes(Phenotype, Cells_colony, color = Phenotype),
  #            size = 5, alpha = 0.5) +
  geom_errorbar(data = extinction_T_av,
                mapping = aes(x = Phenotype, ymax = Cells_colony_mean + Cells_colony_se,
                              ymin = Cells_colony_mean - Cells_colony_se,
                              color = Phenotype),
                width = 0.2) +
  geom_point(data = extinction_T_av,
                mapping = aes(Phenotype, Cells_colony_mean, color = Phenotype, shape = Phenotype),
                size = exp_point_size) +
  scale_color_manual(values = c("Anc" = p_Anc,
                                "Mut" = p_Mut)) +
  scale_shape_manual(values = c("Anc" = sh_Anc,
                                "Mut" = sh_Mut)) +
  scale_y_continuous(trans = "log10", name = "Colony size (CFU)",
                     limits = c(5e6, 5e8),
                     breaks = c(1e7,1e8),
                     labels = sapply(c(7,8),function(i){parse(text = sprintf("10^%d",i))})) +
  fig_aes +
  theme(axis.title.x = element_blank())

### Extinction plot ----
e3 <- ggplot() +
  # geom_point(data = extinction_T, 
  #            mapping = aes(Phenotype, Extinction, color = Phenotype),
  #            size = 5, alpha = 0.5) +
  geom_errorbar(data = extinction_T_av,
                mapping = aes(x = Phenotype, ymax = Extinction_mean + Extinction_se,
                              ymin = Extinction_mean - Extinction_se,
                              color = Phenotype),
                width = 0.2) +
  geom_point(data = extinction_T_av,
             mapping = aes(Phenotype, Extinction_mean, color = Phenotype, shape = Phenotype),
             size = exp_point_size) +
  scale_color_manual(values = c("Anc" = p_Anc,
                                "Mut" = p_Mut)) +
  scale_shape_manual(values = c("Anc" = sh_Anc,
                                "Mut" = sh_Mut)) +
  fig_aes +
  theme(axis.title.x = element_blank())

# D&R Extinction calculation ----
extinction_DR <- plating_DR %>%
  select(Experiment, Phenotype, Replicate, Plate_type, Count) %>%
  mutate(Plate_type = ifelse(Plate_type == "LB", "LB", "Selective")) %>%
  pivot_wider(names_from = Plate_type, values_from = Count) %>%
  mutate(Extinction = 1 - (Selective/LB)) %>%
  # Replace negative values with 0
  mutate(Extinction = ifelse(Extinction < 0, 0, Extinction)) %>%
  select(-LB, -Selective) %>%
  mutate(Host = ifelse(Phenotype == "F", "nalR", "rifR")) %>%
  mutate(Plate_type = ifelse(Phenotype == "F", "Nal50", "Tet50/Rif75"))

# Average values
extinction_DR_av <- extinction_DR %>%
  group_by(Host, Phenotype, Plate_type) %>%
  summarise(Extinction_mean = mean(Extinction),
            Extinction_sd = sd(Extinction),
            n = n()) %>%
  ungroup() %>%
  mutate(Extinction_se = Extinction_sd/sqrt(n))

### Plot DR Extinctions ----
e4 <- ggplot() +
  geom_point(data = extinction_DR, 
             mapping = aes(Phenotype, Extinction, color = Phenotype, shape = Phenotype),
             size = exp_point_size, alpha = 0.4) +
  geom_errorbar(data = extinction_DR_av,
                mapping = aes(x = Phenotype, ymax = Extinction_mean + Extinction_se,
                              ymin = Extinction_mean - Extinction_se,
                              color = Phenotype),
                width = 0.2) +
  geom_point(data = extinction_DR_av,
             mapping = aes(Phenotype, Extinction_mean, color = Phenotype, shape = Phenotype),
             size = exp_point_size) +
  scale_color_manual(values = c("Anc" = p_Anc,
                                "Mut" = p_Mut,
                                "F" = p_F)) +
  scale_shape_manual(values = c("Anc" = sh_Anc,
                                "Mut" = sh_Mut,
                                "F" = sh_F)) +
  fig_aes +
  theme(axis.title.x = element_blank())

# Combine datasets ----
extinction_out <- extinction_T %>%
  bind_rows(extinction_DR)

extinction_out_av <- extinction_T_av %>%
  bind_rows(extinction_DR_av)

### Plot all extinctions together, by media ----
e5 <- ggplot() +
  geom_jitter(data = extinction_out, 
             mapping = aes(Plate_type, Extinction, color = Phenotype, shape = Phenotype),
             size = exp_point_size, alpha = 0.5, width = 0.05) +
  geom_errorbar(data = extinction_out_av,
                mapping = aes(x = Plate_type, ymax = Extinction_mean + Extinction_se,
                              ymin = Extinction_mean - Extinction_se,
                              color = Phenotype),
                width = 0.1) +
  geom_point(data = extinction_out_av,
             mapping = aes(Plate_type, Extinction_mean, color = Phenotype, shape = Phenotype),
             size = exp_point_size) +
  scale_y_continuous(limits = c(-0.05,1)) +
  xlab("Media") +
  scale_color_manual(values = c("Anc" = p_Anc,
                                "Mut" = p_Mut,
                                "F" = p_F)) +
  scale_shape_manual(values = c("Anc" = sh_Anc,
                                "Mut" = sh_Mut,
                                "F" = sh_F)) +
  fig_aes


# Export processed files and plots ----
# csvs
write_csv(extinction_out, paste(output_dir, "extinction_size.csv", sep = "/"))
write_csv(extinction_out_av, paste(output_dir, "extinction_size_av.csv", sep = "/"))

# plots
ggsave(paste(output_dir, "colony_size_T_out.pdf", sep = "/"), s3,
       width = 3, height = 2.5, units = "in")
ggsave(paste(output_dir, "extinction_T_out.pdf", sep = "/"), e3,
       width = 3, height = 2.5, units = "in")

ggsave(paste(output_dir, "extinction_DR_out.pdf", sep = "/"), e4,
       width = 3, height = 2.5, units = "in")
ggsave(paste(output_dir, "extinction_all_out.pdf", sep = "/"), e5,
       width = 5, height = 3, units = "in")

# R objects
save(s3, file = paste(output_dir, "colony_size_T_out.rdata", sep = "/"))
save(e3, file = paste(output_dir, "extinction_T_out.rdata", sep = "/"))
save(e5, file = paste(output_dir, "extinction_all_out.rdata", sep = "/"))
