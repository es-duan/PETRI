# Process and plot Extinction and Colony Size data

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

## Plot all points by colony size ----
### Size plot ----
s1 <- ggplot() +
  geom_jitter(data = processed, aes(Sizes, Cells_colony, color = Treatment),
              alpha = 0.5, width = 0.1, size = 5) +
  scale_y_continuous(trans = "log10", name = "Colony size (CFU)",
                     breaks = c(1e6,1e7,1e8),
                     labels = sapply(c(6,7,8),function(i){parse(text = sprintf("10^%d",i))})) +
  fig_aes

### Extinction plot ----
e1 <- ggplot() +
  geom_jitter(data = processed, aes(Sizes, Extinction, color = Treatment),
              alpha = 0.5, width = 0.1, size = 5) +
  fig_aes

## Compare averages for three groupings ----
# all_points: all points by initial size
# monocultures: monoculture data only
# size_maintained: points where size was maintained
all <- processed %>%
  mutate(group = "all_points")
mono <- processed %>%
  filter(Treatment %in% c("Anc-mono", "Mut-mono")) %>%
  mutate(group = "monocultures")
maintain <- processed %>%
  filter(Sizes %in% c("small-small", "large-large")) %>%
  mutate(group = "size_maintained")
compare <- rbind(all, mono, maintain)

### Size plot ----
s2 <- ggplot(data = compare, aes(Phenotype, Cells_colony, color = Phenotype)) +
  geom_point(alpha = 0.5, size = 5) +
  stat_summary(fun.data="mean_se",geom="errorbar", width=0.2) +
  stat_summary(fun="mean",geom="point", shape = 18, size=5) +
  scale_y_continuous(trans = "log10", name = "Colony size (CFU)",
                     breaks = c(1e6,1e7,1e8),
                     labels = sapply(c(6,7,8),function(i){parse(text = sprintf("10^%d",i))})) +
  facet_grid(~group) +
  scale_color_manual(values = c("Anc" = p_Anc,
                                "Mut" = p_Mut)) +
  fig_aes

### Extinction plot ----
e2 <- ggplot(data = compare, aes(Phenotype, Extinction, color = Phenotype)) +
  geom_point(alpha = 0.5, size = 5) +
  stat_summary(fun.data="mean_se",geom="errorbar", width=0.2) +
  stat_summary(fun="mean",geom="point", shape = 18, size=5) +
  facet_grid(~group) +
  scale_color_manual(values = c("Anc" = p_Anc,
                                "Mut" = p_Mut)) +
  fig_aes

# Move forward with colonies that maintained their sizes ----
extinction_T <- maintain %>%
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
  geom_point(data = extinction_T, 
             mapping = aes(Phenotype, Cells_colony, color = Phenotype),
             size = 5, alpha = 0.5) +
  geom_errorbar(data = extinction_T_av,
                mapping = aes(x = Phenotype, ymax = Cells_colony_mean + Cells_colony_se,
                              ymin = Cells_colony_mean - Cells_colony_se,
                              color = Phenotype),
                width = 0.2) +
  geom_point(data = extinction_T_av,
                mapping = aes(Phenotype, Cells_colony_mean, color = Phenotype),
                size = 5, shape = 18) +
  scale_color_manual(values = c("Anc" = p_Anc,
                                "Mut" = p_Mut)) +
  scale_y_continuous(trans = "log10", name = "Colony size (CFU)",
                     breaks = c(1e6,1e7,1e8),
                     labels = sapply(c(6,7,8),function(i){parse(text = sprintf("10^%d",i))})) +
  fig_aes +
  theme(axis.title.x = element_blank())

### Extinction plot ----
e3 <- ggplot() +
  geom_point(data = extinction_T, 
             mapping = aes(Phenotype, Extinction, color = Phenotype),
             size = 5, alpha = 0.5) +
  geom_errorbar(data = extinction_T_av,
                mapping = aes(x = Phenotype, ymax = Extinction_mean + Extinction_se,
                              ymin = Extinction_mean - Extinction_se,
                              color = Phenotype),
                width = 0.2) +
  geom_point(data = extinction_T_av,
             mapping = aes(Phenotype, Extinction_mean, color = Phenotype),
             size = 5, shape = 18) +
  scale_color_manual(values = c("Anc" = p_Anc,
                                "Mut" = p_Mut)) +
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
             mapping = aes(Phenotype, Extinction, color = Phenotype),
             size = 5, alpha = 0.5) +
  geom_errorbar(data = extinction_DR_av,
                mapping = aes(x = Phenotype, ymax = Extinction_mean + Extinction_se,
                              ymin = Extinction_mean - Extinction_se,
                              color = Phenotype),
                width = 0.2) +
  geom_point(data = extinction_DR_av,
             mapping = aes(Phenotype, Extinction_mean, color = Phenotype),
             size = 5, shape = 18) +
  scale_color_manual(values = c("Anc" = p_Anc,
                                "Mut" = p_Mut,
                                "F" = p_F)) +
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
             mapping = aes(Plate_type, Extinction, color = Phenotype),
             size = 5, alpha = 0.5, width = 0.05) +
  geom_errorbar(data = extinction_out_av,
                mapping = aes(x = Plate_type, ymax = Extinction_mean + Extinction_se,
                              ymin = Extinction_mean - Extinction_se,
                              color = Phenotype),
                width = 0.1) +
  geom_point(data = extinction_out_av,
             mapping = aes(Plate_type, Extinction_mean, color = Phenotype),
             size = 5, shape = 18) +
  scale_y_continuous(limits = c(-0.05,1)) +
  xlab("Media") +
  scale_color_manual(values = c("Anc" = p_Anc,
                                "Mut" = p_Mut,
                                "F" = p_F)) +
  fig_aes


# Export processed files and plots ----
# csvs
write_csv(extinction_out, paste(output_dir, "extinction_size.csv", sep = "/"))
write_csv(extinction_out_av, paste(output_dir, "extinction_size_av.csv", sep = "/"))

# plots
ggsave(paste(output_dir, "colony_size_T_comp.pdf", sep = "/"), s2,
       width = 12, height = 5, units = "in")
ggsave(paste(output_dir, "extinction_T_comp.pdf", sep = "/"), e2,
       width = 12, height = 5, units = "in")

ggsave(paste(output_dir, "colony_size_T_out.pdf", sep = "/"), s3,
       width = 5, height = 4, units = "in")
ggsave(paste(output_dir, "extinction_T_out.pdf", sep = "/"), e3,
       width = 5, height = 4, units = "in")

ggsave(paste(output_dir, "extinction_DR_out.pdf", sep = "/"), e4,
       width = 5, height = 4, units = "in")
ggsave(paste(output_dir, "extinction_all_out.pdf", sep = "/"), e5,
       width = 9, height = 4.5, units = "in")

# R objects
save(s3, file = paste(output_dir, "colony_size_T_out.rdata", sep = "/"))
save(e3, file = paste(output_dir, "extinction_T_out.rdata", sep = "/"))
save(e5, file = paste(output_dir, "extinction_all_out.rdata", sep = "/"))
