# Reformat simrum case study outputs for plotting

# Load packages ----
library(tidyverse)
library(argparse)

# Set arguments parser inputs ----
parser <- ArgumentParser()
parser$add_argument("-t","--treatment", type = "character", help = "Specify Treatment ID")
args <- parser$parse_args()

treatment <- args$treatment

# Load simrun files ----
treatment_folder <- paste("results", "case_study_sims", treatment, sep = "/")

sim_results <- read_csv(paste0(treatment_folder, "/", treatment, "_data.csv"))
sim_results_long <- read_csv(paste0(treatment_folder, "/", treatment, "_data_long.csv"))
sim_phases <- read_csv(paste0(treatment_folder, "/", treatment, "_phases.csv"))

# Reformat data for plotting ----
## Density plot ----
sim_dens <- sim_results_long %>%
  # Remove nutrient values
  filter(Host != "C") %>%
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
  filter(Time %% 1 == 0)

## Phase colors ----
phases_rect <- sim_phases %>%
  group_by(Cycle, Phase) %>%
  summarise(T_start = min(Time),
            T_end = max(Time)) %>%
  ungroup() %>%
  arrange(T_start) %>%
  mutate(Phase = case_when(Phase == "conj" ~ "Conjugation",
                           Phase == "growth" ~ "Growth",
                           Phase == "growout" ~ "Growout",
                           Phase == "tselect" ~ "Transconjugant selection",
                           Phase == "immigration" ~ "Immigration",
                           Phase == "selection" ~ "Plasmid selection"))

## Ratio plots ----

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

# Do not change host identities for protocols without a host switch
if(str_detect(treatment, "E") == FALSE){
  sim_phases2 <- sim_phases %>%
    mutate(Anc = assign_anc(Cycle, Phase, Type, A1, A2)) %>%
    mutate(Mut = assign_mut(Cycle, Phase, Type, M1, M2))
} else{
  sim_phases2 <- sim_phases %>%
    mutate(Anc = A1 + A2) %>%
    mutate(Mut = M1 + M2)
}

sim_phases3 <- sim_phases2 %>%
  mutate(A_M = Anc/Mut,
         M_A = Mut/Anc) %>%
  mutate(Total_plasmid = Anc + Mut) %>%
  mutate(Por_A = Anc/Total_plasmid,
         Por_M = Mut/Total_plasmid)

sim_freq <- sim_phases3 %>%
  select(Cycle, Phase, Time, Por_A, Por_M) %>%
  pivot_longer(cols = -c(Cycle, Phase, Time),
               names_to = "Genotype",
               values_to = "Frequency") %>%
  mutate(Genotype = ifelse(Genotype=="Por_A", "Ancestor", "Mutant"))

# Save files ----
write_csv(sim_dens, paste0(treatment_folder, "/", treatment, "_density_plot_df.csv"))
write_csv(sim_freq, paste0(treatment_folder, "/", treatment, "_frequency_plot_df.csv"))
write_csv(phases_rect, paste0(treatment_folder, "/", treatment, "_phases_plot_df.csv"))
