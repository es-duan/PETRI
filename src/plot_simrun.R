# Produce simrun case study plots
# snakemake ver

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

# Get treatment
treatment <- args$treatment

# Read in files ----
treatment_folder <- paste("results", "case_study_sims", treatment, sep = "/")

sim_dens <- read_csv(paste0(treatment_folder, "/", treatment, "_density_plot_df.csv"))
sim_freq <- read_csv(paste0(treatment_folder, "/", treatment, "_frequency_plot_df.csv"))
phases_rect <- read_csv(paste0(treatment_folder, "/", treatment, "_phases_plot_df.csv")) %>%
  mutate(Phase = factor(Phase, levels = c("Growth", "Conjugation", "Transconjugant selection", "Immigration")))

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

# Plot densities over time ----
## Split values to show strains that are being selected against ----
if(str_detect(treatment, "Dim") == TRUE){
  sim_densL <- sim_dens[1,]
  sim_densB <- sim_dens %>% filter(Density != 0)
  sim_densP <- sim_dens[1,]
  point_size <- 0.01
  
} else{
  sim_dens2 <- sim_dens %>%
    filter(Density != 0) %>%
    mutate(selection = case_when(Phase == "conj" & Cycle %in% c(1,3) & Cell_types %in% c("A1","M1","F2") ~ "yes",
                                 Phase == "conj" & Cycle == 2 & Cell_types %in% c("A2","M2","F1") ~ "yes",
                                 TRUE ~ "no"))
  
  sim_densL <- sim_dens2 %>% filter(selection == "yes")
  sim_densP <- sim_densL %>%
    group_by(Cycle, Cell_types) %>%
    slice_tail(n = 1) %>%
    ungroup()
  point_size <- 2
  
  if(str_detect(treatment, "HFC") == TRUE){
    sim_densB <- filter(sim_dens2, selection == "no")
  }else{
    sim_densB <- sim_densL %>%
      group_by(Cycle, Cell_types) %>%
      slice_head(n = 1) %>%
      ungroup() %>%
      mutate(Time = Time + 0.001) %>%
      rbind(filter(sim_dens2, selection == "no")) 
  }
}


## Plot ----
p1 <- ggplot() + 
  geom_rect(data = phases_rect,
            mapping = aes(xmin = T_start, xmax = T_end, ymin = 1, ymax = Inf, fill = Phase)) +
  scale_fill_manual(values = c("Growth" = p_growth,
                               "Conjugation" = p_conj,
                               "Transconjugant selection" = p_tselect,
                               "Immigration" = p_imm)) +
  geom_line(data = sim_densL,
            mapping = aes(x = Time, y = Density,
                          color = Genotype_plot, linetype = Host,
                          group = interaction(Cycle, Genotype_plot, Host)),
            linewidth = 2, alpha = 0.5) +
  geom_point(data = sim_densP,
             mapping = aes(x = Time, y = Density, color = Genotype_plot),
             shape = 21, size = point_size, fill = "white", alpha = 0.75) +
  geom_line(data = sim_densB,
            mapping = aes(x = Time, y = Density,
                          color = Genotype_plot, linetype = Host,
                          group = interaction(Cycle, Genotype_plot, Host)),
            linewidth = 2) +
  scale_color_manual(values = c("Ancestor" = p_Anc,
                                "Mutant" = p_Mut,
                                "Plasmid-free" = p_F),
                     breaks = c("Ancestor", "Mutant", "Plasmid-free"),
                     labels = c("Ancestor", "Mutant", "Plasmid-free")) +
  scale_linetype_manual(values = c("rifR" = rifR_l,
                                   "nalR" = nalR_l)) +
  geom_hline(yintercept = 1.5, color = "white") +
  labs(x = "Time (hr)",
       y = "Density (CFU/mL)",
       color = "Genotype",
       fill = "Phase",
       linetype = "Host") +
  scale_y_continuous(trans = "log10",
                     breaks=c(1e0,1e2,1e4,1e6,1e8),
                     labels=sapply(c(0,2,4,6,8),function(i){parse(text = sprintf("10^%d",i))}),
                     expand = c(0.01, 0.01),
                     limits = c(0.9, 1.5e9)) +
  scale_x_continuous(expand = c(0.001, 0.001)) +
  guides(color = guide_legend(order=1, override.aes = list(shape = NA)),
         fill = guide_legend(order=2),
         linetype = guide_legend(order=3)) +
  fig_aes

# Save plot
ggsave(paste0(treatment_folder, "/", treatment, "_density_plot.pdf"),
       p1, height = 5, width = 20, units = "in")
saveRDS(p1, paste0(treatment_folder, "/", treatment, "_density_plot.rds"))

# Plot frequencies over time ----
# Use different scales for the plots
if(min(sim_freq$Frequency) < 1e-5){
  f_limits <- c(1e-7, 1)
  f_breaks <- c(1e-6, 1e-4, 1e-2, 1)
  f_labels <- c(-6, -4, -2)
  
} else if(min(sim_freq$Frequency) > 1e-5){
  f_limits <- c(1e-5, 1)
  f_breaks <- c(1e-4, 1e-2, 1)
  f_labels <- c(-4, -2)
}

p2 <- ggplot() + 
  geom_rect(data = phases_rect,
            mapping = aes(xmin = T_start, xmax = T_end, ymin = 0, ymax = Inf, fill = Phase)) +
  scale_fill_manual(values = c("Growth" = p_growth,
                               "Conjugation" = p_conj,
                               "Transconjugant selection" = p_tselect,
                               "Immigration" = p_imm)) +
  geom_line(data = sim_freq,
            mapping = aes(Time, Frequency, color = Genotype),
            linewidth = 2) +
  scale_color_manual(values = c("Ancestor" = p_Anc,
                                "Mutant" = p_Mut)) +
  scale_y_continuous(limits = f_limits,
                     breaks = f_breaks,
                     labels = c(sapply(f_labels,function(i){parse(text = sprintf("10^%d",i))}),1),
                     expand = c(0.01, 0.01),
                     trans = "log10") +
  scale_x_continuous(expand = c(0.001, 0.001)) +
  labs(x = "Time (hr)",
       y = "Frequency",
       color = "Genotype",
       fill = "Phase") +
  guides(color = guide_legend(order=1),
         fill = guide_legend(order=2)) +
  fig_aes

# Save plot
ggsave(paste0(treatment_folder, "/", treatment, "_frequency_plot.pdf"),
       p2, height = 5, width = 20, units = "in")
saveRDS(p2, paste0(treatment_folder, "/", treatment, "_frequency_plot.rds"))
