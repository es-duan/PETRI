# Reformat and plot simrum case study outputs

# Load packages ----
library(tidyverse)
library(colorspace)
library(argparse)
library(cowplot)

# Set arguments parser inputs ----
parser <- ArgumentParser()
parser$add_argument("-t","--treatment", type = "character", help = "Specify Treatment ID")
args <- parser$parse_args()

treatment <- args$treatment

# Set common aesthetics ----
# Colors
p_Anc <- "#8394F6"
p_Mut <- "#8A407A"
p_F <- "gray40"
p_C <- "gray80"
p_batch <- lighten(lighten("#FFEEBD",0.4),0.6)
p_conj <- lighten(lighten("#526AB4",0.4),0.6)
p_T <- lighten(lighten("#E6C5EE",0.15),0.6)

# Plot settings
fig_aes <- theme_cowplot() +
  theme(text = element_text(size=24),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background =element_rect(fill="gray90"),
        strip.text.x = element_text(size = 24),
        axis.title.x = element_text(size=30),
        axis.title.y = element_text(size=30),
        axis.text.x = element_text(size=24),
        axis.text.y = element_text(size=24))

# Load simrun files ----
treatment_folder <- paste("results", "case_study_sims", treatment, sep = "/")

sim_results <- read_csv(paste0(treatment_folder, "/", treatment, "_data.csv"))
sim_results_long <- read_csv(paste0(treatment_folder, "/", treatment, "_data_long.csv"))
sim_phases <- read_csv(paste0(treatment_folder, "/", treatment, "_phases.csv"))

# Reformat data for plotting ----
sim_plot <- sim_results_long %>%
  # Rename values for final plots
  mutate(Host = case_when(Host == "1" ~ "rifR",
                          Host == "2" ~ "nalR",
                          Host == "C" ~ "C")) %>%
  mutate(Genotype = case_when(Genotype == "A" ~ "Ancestor",
                              Genotype == "M" ~ "Mutant",
                              Genotype == "F" ~ "Plasmid-free",
                              Genotype == "C" ~ "C")) %>%
  mutate(Genotype = factor(Genotype, levels = c("Ancestor", "Mutant", "Plasmid-free", "C"))) %>%
  mutate(Genotype_plot = fct_rev(Genotype)) %>%
  mutate(Host = factor(Host, levels = c("rifR", "nalR", "C"))) %>%
  # Decrease number of points for smoother plotting
  group_by(Cycle, Phase, Cell_types, Host, Genotype) %>%
  slice(seq(1, n(), by = 1/0.01*0.02)) %>%
  ungroup()

## Reformat dataset for plotting phase colors ----
phases_rect <- sim_phases %>%
  group_by(Cycle, Phase) %>%
  summarise(T_start = min(Time),
            T_end = max(Time)) %>%
  ungroup() %>%
  arrange(T_start) %>%
  mutate(Phase = case_when(Phase == "conj" ~ "Conjugation",
                           Phase == "growth" ~ "Growth",
                           Phase == "growout" ~ "Growout",
                           Phase == "tselect" ~ "Transconjugant selection"))

p1 <- ggplot() + 
  geom_rect(data = phases_rect,
            mapping = aes(xmin = T_start, xmax = T_end, ymin = 1, ymax = Inf, fill = Phase)) +
  scale_fill_manual(values = c("Growout" = "white",
                               "Growth" = p_batch,
                               "Conjugation" = p_conj,
                               "Transconjugant selection" = p_T)) +
  geom_line(data = sim_plot,
            mapping = aes(x = Time, y = Density,
                          color = Genotype_plot, linetype = Host),
            size = 2) +
  scale_color_manual(values = c("Ancestor" = p_Anc,
                                "Mutant" = p_Mut,
                                "Plasmid-free" = p_F,
                                "C" = p_C),
                     breaks = c("Ancestor", "Mutant", "Plasmid-free", "C"),
                     labels = c("Ancestor", "Mutant", "Plasmid-free", "C")) +
  scale_linetype_manual(values = c("rifR" = "solid",
                                   "nalR" = "dashed",
                                   "C" = "dotted")) +
  geom_hline(yintercept = 1.5, color = "white") +
  labs(x = "Time (hr)",
       y = "Density (CFU/mL)",
       color = "Genotype",
       fill = "Phase",
       linetype = "Host") +
  scale_y_continuous(trans = "log10",limits=c(1,2e9),
                     breaks=c(1e0,1e2,1e4,1e6,1e8),
                     labels=sapply(c(0,2,4,6,8),function(i){parse(text = sprintf("10^%d",i))}),
                     expand = c(0.01, 0.01)) +
  scale_x_continuous(expand = c(0.001, 0.001)) +
  guides(fill = guide_legend(order=1, reverse = TRUE),
         color = guide_legend(order=2),
         linetype = guide_legend(order=3)) +
  fig_aes

# Save plot
p_width <- max(sim_results$Time)/7
p_width <- ifelse(p_width > 20, 20, p_width)
ggsave(paste0(treatment_folder, "/", treatment, "_density_plot.pdf"),
       p1, height = 5, width = p_width, units = "in")


# Plot ratios over time ----

# Functions for designating the value to calculate the ratio with
# For the beginning of a phase, this is the type that was selected for
# For the end of a phase, this is the type that will be selected for
assign_anc <- function(Cycle, Phase, Type, A1, A2){
  anc_out <- c()
  for(i in 1:length(Phase)){
    cycle = Cycle[i]
    phase = Phase[i]
    type = Type[i]
    if (cycle == 0){
      anc <- A1[i]
    } else if (cycle %% 2 == 0){
      # For even cycles, host 2 is the donor and host 1 is the transconjugant
      if (phase == "growth"){
        anc <- A2[i]
      } else if (phase == "conj"){
        anc <- ifelse(type == "start", A2[i], A1[i])
      } else if (phase == "tselect"){
        anc <- A1[i]
      }
      
    } else {
      # For odd cycles, host 1 is the donor and host 2 is the transconjugant
      if (phase == "growth"){
        anc <- A1[i]
      } else if (phase == "conj"){
        anc <- ifelse(type == "start", A1[i], A2[i])
      } else if (phase == "tselect"){
        anc <- A2[i]
      }
    }
    anc_out <- c(anc_out, anc)
  }
  return(anc_out)
}

assign_mut <- function(Cycle, Phase, Type, M1, M2){
  mut_out <- c()
  for(i in 1:length(Phase)){
    cycle = Cycle[i]
    phase = Phase[i]
    type = Type[i]
    if (cycle == 0){
      mut <- M1[i]
    } else if (cycle %% 2 == 0){
      # For even cycles, host 2 is the donor and host 1 is the transconjugant
      if (phase == "growth"){
        mut <- M2[i]
      } else if (phase == "conj"){
        mut <- ifelse(type == "start", M2[i], M1[i])
      } else if (phase == "tselect"){
        mut <- M1[i]
      }
      
    } else {
      # For odd cycles, host 1 is the donor and host 2 is the transconjugant
      if (phase == "growth"){
        mut <- M1[i]
      } else if (phase == "conj"){
        mut <- ifelse(type == "start", M1[i], M2[i])
      } else if (phase == "tselect"){
        mut <- M2[i]
      }
    }
    mut_out <- c(mut_out, mut)
  }
  return(mut_out)
}

sim_plot2 <- sim_phases %>%
  mutate(Anc = assign_anc(Cycle, Phase, Type, A1, A2)) %>%
  mutate(Mut = assign_mut(Cycle, Phase, Type, M1, M2)) %>%
  mutate(A_M = Anc/Mut,
         M_A = Mut/Anc) %>%
  select(Cycle, Phase, Time, A_M, M_A) %>%
  pivot_longer(cols = -c(Cycle, Phase, Time),
               names_to = "Ratio_type",
               values_to = "Ratio") %>%
  mutate(Ratio_type = ifelse(Ratio_type=="A_M", "Ancestor:Mutant", "Mutant:Ancestor"))

# Set ranges for y-axis breaks
p2 <- ggplot() + 
  geom_rect(data = phases_rect,
             mapping = aes(xmin = T_start, xmax = T_end, ymin = 0, ymax = Inf, fill = Phase)) +
  scale_fill_manual(values = c("Growout" = "white",
                               "Growth" = p_batch,
                               "Conjugation" = p_conj,
                               "Transconjugant selection" = p_T)) +
  geom_line(data = sim_plot2,
            aes(Time, Ratio, color = fct_rev(Ratio_type)),
            size = 2) +
  scale_color_manual(values = c("Ancestor:Mutant" = p_Anc,
                                "Mutant:Ancestor" = p_Mut),
                     name = "Ratio",
                     breaks = c("Ancestor:Mutant", "Mutant:Ancestor"),
                     labels = c("Ancestor:Mutant", "Mutant:Ancestor")) +
  scale_y_continuous(trans = "log10") +
  fig_aes

# Save plot
ggsave(paste0(treatment_folder, "/", treatment, "_ratio_plot.pdf"),
       p2, height = 5, width = p_width, units = "in")

# Save plotting files ----
write_csv(sim_plot, paste0(treatment_folder, "/", treatment, "_density_plot_df.csv"))
write_csv(sim_plot2, paste0(treatment_folder, "/", treatment, "_ratio_plot_df.csv"))