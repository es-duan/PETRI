# Script to run invasion experiments of the case study model

# Load packages ----
library(deSolve)
library(tidyverse)
library(cowplot)
library(argparse)
library(jsonlite)

# Set arguments parser inputs ----
parser <- ArgumentParser()
parser$add_argument("-t","--treatment", type = "character", help = "Specify Treatment ID")
parser$add_argument("-c","--colors", help = "JSON string of plot colors")

# Parse arguments
args <- parser$parse_args()

# Get treatment
treatment <- args$treatment

# Load global variables ----
## Colors ----
plot_colors <- fromJSON(args$colors)
p_Anc <- plot_colors[["p_Anc"]]
p_Mut <- plot_colors[["p_Mut"]]
p_F <- plot_colors[["p_F"]]

# Load treatment file ----
treatment_csv <- read.csv("input_data/Treatments_case_study.csv", header = F)

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
Volume = as.numeric(treatment_csv$Volume[row_number])
Dilution_cutoff = 1/Volume
T_volume = as.numeric(treatment_csv$T_volume[row_number])

Protocol = str_split_1(toString(treatment_csv$Protocol[row_number]), pattern = ",")
Protocol_phases = as.numeric(str_split_1(toString(treatment_csv$Protocol_phases[row_number]), pattern = ","))
Selection_type = as.character(treatment_csv$Selection_type[row_number])

Hours_growth = as.numeric(treatment_csv$Hours_growth[row_number])
Hours_conjugation = as.numeric(treatment_csv$Hours_conjugation[row_number])
Hours_t_selection = as.numeric(treatment_csv$Hours_t_selection[row_number])

A_colony = as.numeric(treatment_csv$A_colony[row_number])
M_colony = as.numeric(treatment_csv$M_colony[row_number])

Final_density_t_selection = as.numeric(treatment_csv$Final_density_t_selection[row_number])

Dilution_growth = as.numeric(treatment_csv$Dilution_growth[row_number])
Dilution_conjugation = as.numeric(treatment_csv$Dilution_conjugation[row_number])
Dilution_t_selection = as.numeric(treatment_csv$Dilution_t_selection[row_number])

A1_migrants = as.numeric(treatment_csv$A1_migrants[row_number])
M1_migrants = as.numeric(treatment_csv$M1_migrants[row_number])
F1_migrants = as.numeric(treatment_csv$F1_migrants[row_number])
A2_migrants = as.numeric(treatment_csv$A2_migrants[row_number])
M2_migrants = as.numeric(treatment_csv$M2_migrants[row_number])
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
# Dataframe for saving start and end points of phases
phases <- data.frame()

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
    mutate(Cycle = 0) %>%
    mutate(Phase = "growout") 
  
  # Save start and end times of the phase
  df_start <- head(results, 1) %>%
    mutate(Type = "start")
  df_end <- tail(results, 1) %>%
    mutate(Type = "end")
  
  phases <- rbind(phases,
                  df_start,
                  df_end)
  
} else if(Growout == "no"){
  # Start with initial densities
  results <- data.frame(time = 0, A1 = A1_0, M1 = M1_0, F1 = F1_0,
                        A2 = A2_0, M2 = M2_0, F2 = F2_0, C = C_0,
                        Cycle = NA, Phase = NA)
}

## Define protocol ----
protocol <- c()
for(r in 1:length(Protocol)){
  phase <- rep(Protocol[r], Protocol_phases[r])
  protocol <- c(protocol, phase)
}

## Set densities for each phase ----
for (c in 1:Cycles){
  cycle_df <- tail(results, 1) %>%
    select(-Cycle)
  
  c_phases <- data.frame()
  for (p in 1:length(protocol)){
    # Set phase of protocol
    phase = protocol[p]
    last_phase = ifelse(p == 1, protocol[length(protocol)], protocol[p-1])
    
    if (phase == "growth"){
      ### Growth phase: transfer in plasmid-selecting antibiotics ----
      # Set time period
      time_start = tail(cycle_df$time, 1)
      time_end = time_start + Hours_growth
      times = seq(time_start, time_end, by = Timestep)
      
      # Set densities for the start of the phase
      A1 = tail(cycle_df$A1, 1)
      M1 = tail(cycle_df$M1, 1)
      # Plasmid-free cells are killed
      F1 = 0
      A2 = tail(cycle_df$A2, 1)
      M2 = tail(cycle_df$M2, 1)
      F2 = 0
      
      if (c == 1 & p == 1){
        # Do not dilute strains if this is the first phase of the protocol
        
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
      
      # Save start and end times of the phase
      df_start <- head(out_df, 1) %>%
        mutate(Type = "start")
      df_end <- tail(out_df, 1) %>%
        mutate(Type = "end")
      
      c_phases <- rbind(c_phases,
                        df_start,
                        df_end)
      
      # Add output to cycle dataframe
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
      
      # Save start and end times of the phase
      df_start <- head(out_df, 1) %>%
        mutate(Type = "start")
      df_end <- tail(out_df, 1) %>%
        mutate(Type = "end")
      
      c_phases <- rbind(c_phases,
                        df_start,
                        df_end)
      
      # Add output to cycle dataframe
      cycle_df <- rbind(cycle_df, out_df)
      
    } else if (phase == "transconjugant_selection"){
      ### Transconjugant selection phase ----
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
      
      #### Liquid selection: simulate until target density is reached ----
      if(Selection_type == "liquid"){
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
        
        # Save start and end times of the phase
        df_start <- head(t_out, 1) %>%
          mutate(Type = "start")
        df_end <- tail(t_out, 1) %>%
          mutate(Type = "end")
        
        c_phases <- rbind(c_phases,
                          df_start,
                          df_end)
        
        # Add output to cycle dataframe
        cycle_df <- rbind(cycle_df, t_out)
      } else if (Selection_type == "solid"){
        #### Solid selection: Plate a certain volume of the culture on selective media, then estimate colony size from each cell ----
        # Start with the number of cells in the media (colony number)
        time_start = tail(cycle_df, 1)$time
        
        A1_start = ifelse(A1 * T_volume < 1, 0, A1 * T_volume)
        M1_start = ifelse(M1 * T_volume < 1, 0, M1 * T_volume)
        A2_start = ifelse(A2 * T_volume < 1, 0, A2 * T_volume)
        M2_start = ifelse(M2 * T_volume < 1, 0, M2 * T_volume)
        
        # End with the number of cells in all the colonies
        time_end = time_start + Hours_t_selection
        
        A1_end = (A1_start*A_colony)/Volume
        M1_end = (M1_start*M_colony)/Volume
        A2_end = (A2_start*A_colony)/Volume
        M2_end = (M2_start*M_colony)/Volume
        
        # Determine the dilution factor to reach the target density
        total_density = A1_end + M1_end + A2_end + M2_end
        
        dil_factor = total_density/Final_density_t_selection
        
        # Perform dilution
        time_end_dil = time_end
        
        A1_end_dil = A1_end/dil_factor
        M1_end_dil = M1_end/dil_factor
        A2_end_dil = A2_end/dil_factor
        M2_end_dil = M2_end/dil_factor
        
        # Save the three points as the phase data frame
        # t_out <- data.frame("time" = c(time_start, time_end, time_end_dil),
        #                     "A1" = c(A1_start, A1_end, A1_end_dil),
        #                     "M1" = c(M1_start, M1_end, M1_end_dil),
        #                     "F1" = c(0,0,0),
        #                     "A2" = c(A2_start, A2_end, A2_end_dil),
        #                     "M2" = c(M2_start, M2_end, M2_end_dil),
        #                     "F2" = c(0,0,0),
        #                     "C" = c(C_0,C_0,C_0)) %>%
        #   mutate(Phase = "tselect")
        
        # Save two data points (start and diluted end only)
        t_out <- data.frame("time" = c(time_start, time_end_dil),
                            "A1" = c(A1_start, A1_end_dil),
                            "M1" = c(M1_start, M1_end_dil),
                            "F1" = c(0,0),
                            "A2" = c(A2_start, A2_end_dil),
                            "M2" = c(M2_start, M2_end_dil),
                            "F2" = c(0,0),
                            "C" = c(C_0,C_0)) %>%
          mutate(Phase = "tselect")
        
        # Save start and end times of the phase
        df_start <- head(t_out, 1) %>%
          mutate(Type = "start")
        df_end <- tail(t_out, 1) %>%
          mutate(Type = "end")
        
        c_phases <- rbind(c_phases,
                          df_start,
                          df_end)
        
        # Add output to cycle dataframe
        cycle_df <- rbind(cycle_df, t_out)
      }
      
    } else if (phase == "immigration"){
      ### Immigration phase: add migrants ----
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
      
      if (c == 1 & p == 1){
        # Do not dilute strains if this is the first phase of the protocol
        
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
      
      # Add migrants
      A1 = A1 + A1_migrants
      M1 = M1 + M1_migrants
      F1 = F1 + F1_migrants
      A2 = A2 + A2_migrants
      M2 = M2 + M2_migrants
      F2 = F2 + F2_migrants
      
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
        mutate(Phase = "immigration")
      #out_df <- out_df %>% slice(seq(1, n(), by = 1/Timestep*0.01))
      
      # Save start and end times of the phase
      df_start <- head(out_df, 1) %>%
        mutate(Type = "start")
      df_end <- tail(out_df, 1) %>%
        mutate(Type = "end")
      
      c_phases <- rbind(c_phases,
                        df_start,
                        df_end)
      
      # Add output to cycle dataframe
      cycle_df <- rbind(cycle_df, out_df)
      
    } else if (phase == "selection"){
      ### Selection phase: kill plasmid-free cells ----
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
      
      if (c == 1 & p == 1){
        # Do not dilute strains if this is the first phase of the protocol
        
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
      
      # Set plasmid-free cells to 0
      F1 = 0
      F2 = 0
      
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
        mutate(Phase = "selection")
      #out_df <- out_df %>% slice(seq(1, n(), by = 1/Timestep*0.01))
      
      # Save start and end times of the phase
      df_start <- head(out_df, 1) %>%
        mutate(Type = "start")
      df_end <- tail(out_df, 1) %>%
        mutate(Type = "end")
      
      c_phases <- rbind(c_phases,
                        df_start,
                        df_end)
      
      # Add output to cycle dataframe
      cycle_df <- rbind(cycle_df, out_df)
    }
  }
  cycle_df <- cycle_df %>%
    mutate(Cycle = c)
  c_phases <- c_phases %>%
    mutate(Cycle = c)
  
  # Save results to main dataframes
  results <- rbind(results, cycle_df)
  phases <- rbind(phases, c_phases)
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

gg <- ggplot(data = plot_out, aes(x = Time, y = Density, group = interaction(Host, Genotype), color = Genotype)) + 
  geom_line(aes(linetype = Host), size = 1) +
  scale_color_manual(values = c("A" = p_Anc,
                                "M" = p_Mut,
                                "F" = p_F))+
  scale_linetype_manual(values = c("1" = "solid", "2" = "dashed", "C" = "dotted")) +
  scale_y_continuous(trans = "log10",breaks=c(1e1,1e3,1e5,1e7,1e9)) +
  geom_hline(yintercept = 1) +
  theme_bw()
ggsave(paste(treatment_folder, paste(treatment, '_sim_plot.pdf', sep = ""), sep = '/'), gg, width = 20, height = 5, units = "in")


# Save the final dataframes ----
results_out <- results %>%
  drop_na() %>%
  rename(Time = time)

phases_out <- phases %>%
  drop_na() %>%
  rename(Time = time)

write.csv(x = results_out, file = paste(treatment_folder, paste(treatment, '_data.csv', sep = ""), sep = '/'))
write.csv(x = plot_out, file = paste(treatment_folder, paste(treatment, '_data_long.csv', sep = ""), sep = '/'))
write.csv(x = phases_out, file = paste(treatment_folder, paste(treatment, '_phases.csv', sep = ""), sep = '/'))

