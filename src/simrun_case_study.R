# Script to run invasion experiments of the case study model

# Load packages ----
library(deSolve)
library(tidyverse)
library(cowplot)
library(argparse)
theme_set(theme_cowplot())

# Set arguments parser inputs ----
parser <- ArgumentParser()
parser$add_argument("-t","--treatment", type = "character", help = "Specify Treatment ID")
args <- parser$parse_args()

treatment <- args$treatment

# Load treatment file ----
treatment_csv <- read.csv("input_data/Treatments.csv", header = F)

## Transpose data to make it easier for reading in parameters ----
tr_colnames <- treatment_csv[[1]]
treatment_csv <- as.data.frame(t(treatment_csv[,-1]))
colnames(treatment_csv) <- tr_colnames

## Select row with treatment to run ----
row_number <- which(treatment_csv$Treatment_ID == treatment)

## Create folder to store results
output_folder <- "results/case_study_sims"
treatment_folder <- file.path(output_folder, treatment)

if (!dir.exists(treatment_folder)) {
  dir.create(treatment_folder, recursive = TRUE)
} else{
  print(paste("Treatment",treatment,"folder exists. Rewriting previous run."))
}

# Set values for parameters, experiment settings, and initial variables ----
## Parameters ----
Gamma_A1.F1_max  = as.numeric(treatment_csv$gammaA1.F1[row_number])
Gamma_A1.F1_base = as.numeric(treatment_csv$gammaA1.F1_base[row_number])
Gamma_M1.F1_max  = as.numeric(treatment_csv$gammaM1.F1[row_number])
Gamma_M1.F1_base = as.numeric(treatment_csv$gammaM1.F1_base[row_number])
Gamma_A2.F2_max  = as.numeric(treatment_csv$gammaA2.F2[row_number])
Gamma_A2.F2_base = as.numeric(treatment_csv$gammaA2.F2_base[row_number])
Gamma_M2.F2_max  = as.numeric(treatment_csv$gammaM2.F2[row_number])
Gamma_M2.F2_base = as.numeric(treatment_csv$gammaM2.F2_base[row_number])
Gamma_A1.F2_max  = as.numeric(treatment_csv$gammaA1.F2[row_number])
Gamma_A1.F2_base = as.numeric(treatment_csv$gammaA1.F2_base[row_number])
Gamma_M1.F2_max  = as.numeric(treatment_csv$gammaM1.F2[row_number])
Gamma_M1.F2_base = as.numeric(treatment_csv$gammaM1.F2_base[row_number])
Gamma_A2.F1_max  = as.numeric(treatment_csv$gammaA2.F1[row_number])
Gamma_A2.F1_base = as.numeric(treatment_csv$gammaA2.F1_base[row_number])
Gamma_M2.F1_max  = as.numeric(treatment_csv$gammaM2.F1[row_number])
Gamma_M2.F1_base = as.numeric(treatment_csv$gammaM2.F1_base[row_number])
Psi_A1_max = as.numeric(treatment_csv$psiA1[row_number])
Psi_M1_max = as.numeric(treatment_csv$psiM1[row_number])
Psi_F1_max = as.numeric(treatment_csv$psiF1[row_number])
Psi_A2_max = as.numeric(treatment_csv$psiA2[row_number])
Psi_M2_max = as.numeric(treatment_csv$psiM2[row_number])
Psi_F2_max = as.numeric(treatment_csv$psiF2[row_number])
Sigma_A1 = as.numeric(treatment_csv$sigmaA1[row_number])
Sigma_M1 = as.numeric(treatment_csv$sigmaM1[row_number])
Sigma_A2 = as.numeric(treatment_csv$sigmaA2[row_number])
Sigma_M2 = as.numeric(treatment_csv$sigmaM2[row_number])
e_0 = as.numeric(treatment_csv$e[row_number])
Q_0 = as.numeric(treatment_csv$Q[row_number])

## Experimental settings ----
Growout = as.character(treatment_csv$Growout[row_number])
Final_density_growout = as.numeric(treatment_csv$Final_density_growout[row_number])
Timestep = as.numeric(treatment_csv$Timestep[row_number])
Cycles = as.numeric(treatment_csv$Cycles[row_number])
Dilution_cutoff = as.numeric(treatment_csv$Dilution_cutoff[row_number])

Phases_growth = as.numeric(treatment_csv$Phases_growth[row_number])
Phases_conjugation = as.numeric(treatment_csv$Phases_conjugation[row_number])
Phases_t_selection = as.numeric(treatment_csv$Phases_t_selection[row_number])

Hours_growth = as.numeric(treatment_csv$Hours_growth[row_number])
Hours_conjugation = as.numeric(treatment_csv$Hours_conjugation[row_number])

Final_density_t_selection = as.numeric(treatment_csv$Final_density_t_selection[row_number])

Dilution_growth = as.numeric(treatment_csv$Dilution_growth[row_number])
Dilution_conjugation = as.numeric(treatment_csv$Dilution_conjugation[row_number])
Dilution_t_selection = as.numeric(treatment_csv$Dilution_t_selection[row_number])

F1_migrants = as.numeric(treatment_csv$F1_migrants[row_number])
F2_migrants = as.numeric(treatment_csv$F2_migrants[row_number])


## Variables ----
A1_0 = as.numeric(treatment_csv$A1_0[row_number])
M1_0 = as.numeric(treatment_csv$M1_0[row_number])
F1_0 = as.numeric(treatment_csv$F1_0[row_number])
A2_0 = as.numeric(treatment_csv$A2_0[row_number])
M2_0 = as.numeric(treatment_csv$M2_0[row_number])
F2_0 = as.numeric(treatment_csv$F2_0[row_number])
C_0  = as.numeric(treatment_csv$C[row_number])


# Define the model ----
model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    # guarantees that the variables current values can be called
    A1 <- state["A1"]
    M1 <- state["M1"]
    F1 <- state["F1"]
    A2 <- state["A2"]
    M2 <- state["M2"]
    F2 <- state["F2"]
    C  <- state["C"]
    
    # within.species resource.dependent conjugation equations
    gamma_A1.F1_C <- gamma_A1.F1_base + ((gamma_A1.F1_max - gamma_A1.F1_base) * (C / (Q + C)))
    gamma_M1.F1_C <- gamma_M1.F1_base + ((gamma_M1.F1_max - gamma_M1.F1_base) * (C / (Q + C)))
    gamma_A2.F2_C <- gamma_A2.F2_base + ((gamma_A2.F2_max - gamma_A2.F2_base) * (C / (Q + C)))
    gamma_M2.F2_C <- gamma_M2.F2_base + ((gamma_M2.F2_max - gamma_M2.F2_base) * (C / (Q + C)))
    
    # cross.species resource.dependent conjugation equations
    gamma_A1.F2_C <- gamma_A1.F2_base + ((gamma_A1.F2_max - gamma_A1.F2_base) * (C / (Q + C)))
    gamma_M1.F2_C <- gamma_M1.F2_base + ((gamma_M1.F2_max - gamma_M1.F2_base) * (C / (Q + C)))
    gamma_A2.F1_C <- gamma_A2.F1_base + ((gamma_A2.F1_max - gamma_A2.F1_base) * (C / (Q + C)))
    gamma_M2.F1_C <- gamma_M2.F1_base + ((gamma_M2.F1_max - gamma_M2.F1_base) * (C / (Q + C)))
    
    # resource-dependent growth equations
    psi_A1_C <- psi_A1_max * (C / (Q + C))
    psi_M1_C <- psi_M1_max * (C / (Q + C))
    psi_F1_C <- psi_F1_max * (C / (Q + C))    
    psi_A2_C <- psi_A2_max * (C / (Q + C))
    psi_M2_C <- psi_M2_max * (C / (Q + C))
    psi_F2_C <- psi_F2_max * (C / (Q + C))
    
    # system of ordinary differential equations
    dA1 <- (psi_A1_C * (1 - sigma_A1) * A1) + (gamma_A1.F1_C * A1 * F1) + (gamma_A2.F1_C * A2 * F1)
    dM1 <- (psi_M1_C * (1 - sigma_M1) * M1) + (gamma_M1.F1_C * M1 * F1) + (gamma_M2.F1_C * M2 * F1)
    dF1 <- (psi_F1_C * F1) + (psi_A1_C * sigma_A1 * A1) + (psi_M1_C * sigma_M1 * M1) - (gamma_A1.F1_C * A1 * F1) - (gamma_A2.F1_C * A2 * F1) - (gamma_M1.F1_C * M1 * F1) - (gamma_M2.F1_C * M2 * F1)
    dA2 <- (psi_A2_C * (1 - sigma_A2) * A2) + (gamma_A1.F2_C * A1 * F2) + (gamma_A2.F2_C * A2 * F2)
    dM2 <- (psi_M2_C * (1 - sigma_M2) * M2) + (gamma_M1.F2_C * M1 * F2) + (gamma_M2.F2_C * M2 * F2)
    dF2 <- (psi_F2_C * F2) + (psi_A2_C * sigma_A2 * A2) + (psi_M2_C * sigma_M2 * M2) - (gamma_A1.F2_C * A1 * F2) - (gamma_A2.F2_C * A2 * F2) - (gamma_M1.F2_C * M1 * F2) - (gamma_M2.F2_C * M2 * F2)
    dC  <- -(((psi_A1_C * A1) + (psi_M1_C * M1) + (psi_F1_C * F1) + (psi_A2_C * A2) + (psi_M2_C * M2) + (psi_F2_C * F2)) * e)
    
    list(c(dA1, dM1, dF1, dA2, dM2, dF2, dC))
  })
}

## Initial conditions ----
state <- c(A1 = A1_0, 
           M1 = M1_0, 
           F1 = F1_0, 
           A2 = A2_0, 
           M2 = M2_0, 
           F2 = F2_0, 
           C  = C_0)

## Parameters ----
parameters <- c(gamma_A1.F1_max = Gamma_A1.F1_max,
                gamma_A1.F1_base = Gamma_A1.F1_base,
                gamma_M1.F1_max = Gamma_M1.F1_max,
                gamma_M1.F1_base = Gamma_M1.F1_base,
                gamma_A2.F2_max = Gamma_A2.F2_max,
                gamma_A2.F2_base = Gamma_A2.F2_base,
                gamma_M2.F2_max = Gamma_M2.F2_max,
                gamma_M2.F2_base = Gamma_M2.F2_base,
                gamma_A1.F2_max = Gamma_A1.F2_max,
                gamma_A1.F2_base = Gamma_A1.F2_base,
                gamma_M1.F2_max = Gamma_M1.F2_max,
                gamma_M1.F2_base = Gamma_M1.F2_base,
                gamma_A2.F1_max = Gamma_A2.F1_max,
                gamma_A2.F1_base = Gamma_A2.F1_base,
                gamma_M2.F1_max = Gamma_M2.F1_max,
                gamma_M2.F1_base = Gamma_M2.F1_base,
                psi_A1_max = Psi_A1_max,
                psi_M1_max = Psi_M1_max,
                psi_F1_max = Psi_F1_max,  
                psi_A2_max = Psi_A2_max,
                psi_M2_max = Psi_M2_max,
                psi_F2_max = Psi_F2_max, 
                sigma_A1 = Sigma_A1, 
                sigma_M1 = Sigma_M1, 
                sigma_A2 = Sigma_A2, 
                sigma_M2 = Sigma_M2, 
                e = e_0,
                Q = Q_0)


# Simulation ----
## Initial settings ----
if(Growout == "yes"){
  # For full invasion experiments, include a growout phase that assumes the
  # invading type arose during a transfer prior to starting the selection protocol
  
  # Dataframe to store results
  results <- data.frame(time = 0, A1 = A1_0, M1 = M1_0, F1 = F1_0, A2 = A2_0, M2 = M2_0, F2 = F2_0, C = C_0)
  
  init_total_density = A1_0 + M1_0
  while (init_total_density < Final_density_growout){
    # Perform step by step simulation until target density is reached
    time_start = tail(results$time, 1)
    time_end = time_start + Timestep
    times <- c(time_start, time_end)
    
    # Set densities
    A1 = tail(results$A1, 1)
    M1 = tail(results$M1, 1)
    F1 = tail(results$F1, 1)
    A2 = tail(results$A2, 1)
    M2 = tail(results$M2, 1)
    F2 = tail(results$F2, 1)
    C = tail(results$C, 1)
    
    # Euler simulation
    state <- c(A1 = A1, 
               M1 = M1, 
               F1 = F1, 
               A2 = A2, 
               M2 = M2, 
               F2 = F2, 
               C  = C)
    
    out <- ode(y = state, times = times, func = model, parms = parameters, method = "euler", hini = Timestep)
    out_df <- as.data.frame(out)
    #out_df <- out_df %>% slice(seq(1, n(), by = 1/Timestep*0.01))
    results <- rbind(results, out_df)
    
    # Re-calculate total density
    A1 = tail(results$A1, 1)
    M1 = tail(results$M1, 1)
    init_total_density = A1 + M1
    
    # Maintain nutrient levels
    C = tail(results$C, 1)
  }
  results <- unique(results) %>%
    mutate(Phase = "growout")
  
} else if(Growout == "no"){
  # Start with initial densities
  results <- data.frame(time = 0, A1 = A1_0, M1 = M1_0, F1 = F1_0,
                        A2 = A2_0, M2 = M2_0, F2 = F2_0, C = C_0,
                        Cycle = NA, Phase = NA)
  
}

## Define protocol ----
protocol <- c(rep("growth", Phases_growth),
              rep("conjugation", Phases_conjugation),
              rep("transconjugant_selection", Phases_t_selection))

## Set densities for each phase ----
for (c in 1:Cycles){
  cycle_df <- tail(results, 1) %>%
    select(-Cycle)
  for (p in 1:length(protocol)){
    # Set phase of protocol
    phase = protocol[p]
    last_phase = ifelse(p == 1, protocol[length(protocol)], protocol[p-1])
    
    if (phase == "growth"){
      ### Growth phase ----
      # Set time period
      time_start = tail(cycle_df$time, 1)
      time_end = time_start + Hours_growth
      times = seq(time_start, time_end, by = Timestep)
      
      # Set densities for the start of the phase
      A1 = tail(cycle_df$A1, 1)
      M1 = tail(cycle_df$M1, 1)
      F1 = tail(cycle_df$F1, 1)
      A2 = tail(cycle_df$A2, 1)
      M2 = tail(cycle_df$M2, 1)
      F2 = tail(cycle_df$F2, 1)
      
      if (time_start == 0){
        # If simulation is starting with a growth phase, do not dilute strains
        
      } else if (last_phase == "transconjugant_selection"){
        # If transitioning from a transconjugant selection phase, do not dilute strains
        
      } else {
        # Else, dilute strains first
        A1 = ifelse((A1 * Dilution_growth) < Dilution_cutoff, 0, A1 * Dilution_growth)
        M1 = ifelse((M1 * Dilution_growth) < Dilution_cutoff, 0, M1 * Dilution_growth)
        F1 = ifelse((F1 * Dilution_growth) < Dilution_cutoff, 0, F1 * Dilution_growth)
        A2 = ifelse((A2 * Dilution_growth) < Dilution_cutoff, 0, A2 * Dilution_growth)
        M2 = ifelse((M2 * Dilution_growth) < Dilution_cutoff, 0, M2 * Dilution_growth)
        F2 = ifelse((F2 * Dilution_growth) < Dilution_cutoff, 0, F2 * Dilution_growth)
      }
      
      # Euler simulation
      C  = as.numeric(treatment_csv$C[row_number])
      
      state <- c(A1 = A1, 
                 M1 = M1, 
                 F1 = F1, 
                 A2 = A2, 
                 M2 = M2, 
                 F2 = F2, 
                 C  = C)
      
      out <- ode(y = state, times = times, func = model, parms = parameters, method = "euler", hini = Timestep)
      out_df <- as.data.frame(out) %>%
        mutate(Phase = "growth")
      #out_df <- out_df %>% slice(seq(1, n(), by = 1/Timestep*0.01))
      cycle_df <- rbind(cycle_df, out_df)
      
    } else if (phase == "conjugation"){
      ### Conjugation phase ----
      # Set time period
      time_start = tail(cycle_df$time, 1)
      time_end = time_start + Hours_conjugation
      times = seq(time_start, time_end, by = Timestep)
      
      # Set densities for the start of the phase
      A1 = tail(cycle_df$A1, 1)
      M1 = tail(cycle_df$M1, 1)
      F1 = tail(cycle_df$F1, 1)
      A2 = tail(cycle_df$A2, 1)
      M2 = tail(cycle_df$M2, 1)
      F2 = tail(cycle_df$F2, 1)
      
      if (last_phase == "growth"){
        # If entering from a growth phase, perform a dilution
        A1 = ifelse((A1 * Dilution_conjugation) < Dilution_cutoff, 0, A1 * Dilution_conjugation)
        M1 = ifelse((M1 * Dilution_conjugation) < Dilution_cutoff, 0, M1 * Dilution_conjugation)
        F1 = ifelse((F1 * Dilution_conjugation) < Dilution_cutoff, 0, F1 * Dilution_conjugation)
        A2 = ifelse((A2 * Dilution_conjugation) < Dilution_cutoff, 0, A2 * Dilution_conjugation)
        M2 = ifelse((M2 * Dilution_conjugation) < Dilution_cutoff, 0, M2 * Dilution_conjugation)
        F2 = ifelse((F2 * Dilution_conjugation) < Dilution_cutoff, 0, F2 * Dilution_conjugation)
      }
      # Add plasmid-free migrants
      if(c %% 2 == 0){
        # For even cycles, add F1
        F1 = F1_migrants
      } else {
        # For odd cycles, add F2
        F2 = F2_migrants
      }
      
      # Euler simulation
      C  = as.numeric(treatment_csv$C[row_number])
      
      state <- c(A1 = A1, 
                 M1 = M1, 
                 F1 = F1, 
                 A2 = A2, 
                 M2 = M2, 
                 F2 = F2, 
                 C  = C)
      
      out <- ode(y = state, times = times, func = model, parms = parameters, method = "euler", hini = Timestep)
      out_df <- as.data.frame(out) %>%
        mutate(Phase = "conj")
      #out_df <- out_df %>% slice(seq(1, n(), by = 1/Timestep*0.01))
      cycle_df <- rbind(cycle_df, out_df)
      
    } else if (phase == "transconjugant_selection"){
      ### Transconjugant selection phase: simulate until target density is reached ----
      # Set densities for the start of the phase
      if(c %% 2 == 0){
        # For even cycles, select for F1 transconjugants
        A1 = tail(cycle_df$A1, 1)
        M1 = tail(cycle_df$M1, 1)
        F1 = 0
        A2 = 0
        M2 = 0
        F2 = 0
        C  = as.numeric(treatment_csv$C[row_number])
      } else {
        # For odd cycles, select for F2 transconjugants
        A1 = 0
        M1 = 0
        F1 = 0
        A2 = tail(cycle_df$A2, 1)
        M2 = tail(cycle_df$M2, 1)
        F2 = 0
        C  = as.numeric(treatment_csv$C[row_number])
      }
      
      # Perform dilution
      A1 = ifelse((A1 * Dilution_t_selection) < Dilution_cutoff, 0, A1 * Dilution_t_selection)
      M1 = ifelse((M1 * Dilution_t_selection) < Dilution_cutoff, 0, M1 * Dilution_t_selection)
      F1 = ifelse((F1 * Dilution_t_selection) < Dilution_cutoff, 0, F1 * Dilution_t_selection)
      A2 = ifelse((A2 * Dilution_t_selection) < Dilution_cutoff, 0, A2 * Dilution_t_selection)
      M2 = ifelse((M2 * Dilution_t_selection) < Dilution_cutoff, 0, M2 * Dilution_t_selection)
      F2 = ifelse((F2 * Dilution_t_selection) < Dilution_cutoff, 0, F2 * Dilution_t_selection)
      
      # Calculate total density of plasmid-containing types
      total_density = A1 + M1 + A2 + M2
      t_out <- tail(cycle_df, 1) %>%
        select(-Phase)
      while (total_density < Final_density_t_selection){
        # Perform step by step simulation until target density is reached
        time_start = tail(t_out$time, 1)
        time_end = time_start + Timestep
        times <- c(time_start, time_end)
        
        # Euler simulation
        state <- c(A1 = A1, 
                   M1 = M1, 
                   F1 = F1, 
                   A2 = A2, 
                   M2 = M2, 
                   F2 = F2, 
                   C  = C)
        
        out <- ode(y = state, times = times, func = model, parms = parameters, method = "euler", hini = Timestep)
        out_df <- as.data.frame(out)
        #out_df <- out_df %>% slice(seq(1, n(), by = 1/Timestep*0.01))
        t_out <- rbind(t_out, out_df)
        
        # Re-calculate total density of plasmid-containing types
        A1 = tail(t_out$A1, 1)
        M1 = tail(t_out$M1, 1)
        A2 = tail(t_out$A2, 1)
        M2 = tail(t_out$M2, 1)
        total_density = A1 + M1 + A2 + M2
        
        # Maintain nutrient levels
        C = tail(t_out$C, 1)
      }
      # Remove the first row of the dataframe (repeat of results)
      t_out <- t_out[-1,]
      # Remove repeat rows from iterations
      t_out <- unique(t_out) %>%
        mutate(Phase = "tselect")
      # Merge with main dataframe
      cycle_df <- rbind(cycle_df, t_out)
    }
  }
  cycle_df <- cycle_df %>%
    mutate(Cycle = c)
  results <- rbind(results, cycle_df)
  print(paste("Cycle",c,"complete."))
}


# Plot the results ----
plot_out <- results %>% 
  drop_na() %>%
  rename(Time = time) %>%
  pivot_longer(
    cols = -c(Time, Cycle, Phase),  # Select all columns except 'time'
    names_to = "Cell_types",  # New column to hold the former column names
    values_to = "Density") %>% # New column to hold the values
  mutate(Host = ifelse(Cell_types == "C", "C", substr(Cell_types, 2, 2)), 
         Genotype = substr(Cell_types, 1, 1)) %>%
  mutate(Density = ifelse(Density < 1, 0, Density))

write.csv(x = results, file = paste(treatment_folder, paste(treatment, '_data.csv', sep = ""), sep = '/'))
write.csv(x = plot_out, file = paste(treatment_folder, paste(treatment, '_data_long.csv', sep = ""), sep = '/'))

gg <- ggplot(data = plot_out, aes(x = Time, y = Density, group = interaction(Host, Genotype), color = Genotype)) + 
  geom_line(aes(linetype = Host), size = 1) +
  scale_color_manual(values = c("A" = "red", "M" = "blue", "F" = "grey"))+
  scale_linetype_manual(values = c("1" = "solid", "2" = "dashed", "C" = "dotted")) +
  scale_y_continuous(trans = "log10",breaks=c(1e1,1e3,1e5,1e7,1e9)) +
  geom_hline(yintercept = 1)
ggsave(paste(treatment_folder, paste(treatment, '_plot.pdf', sep = ""), sep = '/'), gg, width = 30, height = 7, units = "in")


