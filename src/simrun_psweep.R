# Script to run invasion simulations in parallel for parameter sweep plots

# Load packages ----
library(deSolve)
library(tidyverse)
library(argparse)
library(foreach)
library(doParallel)

# Set arguments parser inputs ----
parser <- ArgumentParser()
parser$add_argument("-p","--psweepsetting", type = "character", help = "Specify parameter sweep setting")
parser$add_argument("-t","--threads", type = "numeric", help = "Specify number of cores given to the script")

# Parse arguments
args <- parser$parse_args()
ps <- args$psweepsetting
n_cores <- args$threads

# Read in files ----
# Specify output folder
output_folder <- paste("results", "parameter_sweeps", ps, sep = "/")

setting_list <- readRDS(paste0(output_folder, "/", ps, "_settings.rds"))
sweep_param <- read_csv(paste0(output_folder, "/", ps, "_params.csv"))

# Set values for experiment settings ----

## Experimental settings ----
Timestep = as.numeric(setting_list$Timestep)
Cycles = as.numeric(setting_list$Cycles)
Volume = as.numeric(setting_list$Volume)
Dilution_cutoff = 1/Volume
T_volume = as.numeric(setting_list$T_volume)
Colony_bottleneck = as.numeric(setting_list$Colony_bottleneck)
Colony_selection = as.character(setting_list$Colony_selection)

Protocol = str_split_1(toString(setting_list$Protocol), pattern = ",")
Protocol_phases = as.numeric(str_split_1(toString(setting_list$Protocol_phases), pattern = ","))
Selection_type = as.character(setting_list$Selection_type)

Hours_growth = as.numeric(setting_list$Hours_growth)
Hours_conjugation = as.numeric(setting_list$Hours_conjugation)
Hours_t_selection = as.numeric(setting_list$Hours_t_selection)

A_colony = as.numeric(setting_list$A_colony)
M_colony = as.numeric(setting_list$M_colony)

Final_density_t_selection = as.numeric(setting_list$Final_density_t_selection)

Dilution_growth = as.numeric(setting_list$Dilution_growth)
Dilution_conjugation = as.numeric(setting_list$Dilution_conjugation)
Dilution_t_selection = as.numeric(setting_list$Dilution_t_selection)

A1_migrants = as.numeric(setting_list$A1_migrants)
M1_migrants = as.numeric(setting_list$M1_migrants)
F1_migrants = as.numeric(setting_list$F1_migrants)
A2_migrants = as.numeric(setting_list$A2_migrants)
M2_migrants = as.numeric(setting_list$M2_migrants)
F2_migrants = as.numeric(setting_list$F2_migrants)
Immigration_ratio = as.numeric(setting_list$Immigration_ratio)


# Define the model ----
model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    # guarantees that the variables current values can be called
    # A1 <- state["A1"]
    # M1 <- state["M1"]
    # F1 <- state["F1"]
    # A2 <- state["A2"]
    # M2 <- state["M2"]
    # F2 <- state["F2"]
    # C  <- state["C"]
    
    # within.species resource.dependent conjugation equations
    gamma_A1.F1_C <- ((gamma_A1.F1_max) * (C / (Q + C)))
    gamma_M1.F1_C <- ((gamma_M1.F1_max) * (C / (Q + C)))
    gamma_A2.F2_C <- ((gamma_A2.F2_max) * (C / (Q + C)))
    gamma_M2.F2_C <- ((gamma_M2.F2_max) * (C / (Q + C)))
    
    # cross.species resource.dependent conjugation equations
    gamma_A1.F2_C <- ((gamma_A1.F2_max) * (C / (Q + C)))
    gamma_M1.F2_C <- ((gamma_M1.F2_max) * (C / (Q + C)))
    gamma_A2.F1_C <- ((gamma_A2.F1_max) * (C / (Q + C)))
    gamma_M2.F1_C <- ((gamma_M2.F1_max) * (C / (Q + C)))
    
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

# 1. Root Function: Detects when C drops to a microscopic threshold
rootfunc <- function(time, state, parameters) {
  return(state["C"] - 1e-10) 
}

# 2. Event Function: Triggers when the root is hit to clamp values cleanly
eventfunc <- function(time, state, parameters) {
  # Clamp C to exactly 0 to stop negative resource drain
  state["C"] <- 0  
  
  # Clamp populations to 0 to catch any tiny floating-point negatives
  state["A1"] <- max(state["A1"], 0)
  state["M1"] <- max(state["M1"], 0)
  state["F1"] <- max(state["F1"], 0)
  state["A2"] <- max(state["A2"], 0)
  state["M2"] <- max(state["M2"], 0)
  state["F2"] <- max(state["F2"], 0)
  
  return(state)
}

# Initiate for loop for parameter sweep ----
## Initiate cluster setup ----

# Register cluster
cluster <- makeCluster(n_cores)
registerDoParallel(cluster)

# For loop ----
sweep_out <- foreach(i = 1:nrow(sweep_param),
        .packages = c("tidyverse", "deSolve"),
        .combine = rbind) %dopar% {
  ## Set parameters ----
  Gamma_A1.F1_max  = as.numeric(setting_list$Ref_gamma)
  Gamma_M1.F1_max  = as.numeric(sweep_param$gamma_M[i])
  Gamma_A2.F2_max  = as.numeric(setting_list$Ref_gamma)
  Gamma_M2.F2_max  = as.numeric(sweep_param$gamma_M[i])
  Gamma_A1.F2_max  = as.numeric(setting_list$Ref_gamma)
  Gamma_M1.F2_max  = as.numeric(sweep_param$gamma_M[i])
  Gamma_A2.F1_max  = as.numeric(setting_list$Ref_gamma)
  Gamma_M2.F1_max  = as.numeric(sweep_param$gamma_M[i])
  Psi_A1_max = as.numeric(setting_list$Ref_psi)
  Psi_M1_max = as.numeric(sweep_param$psi_M[i])
  Psi_F1_max = as.numeric(setting_list$psiF1)
  Psi_A2_max = as.numeric(setting_list$Ref_psi)
  Psi_M2_max = as.numeric(sweep_param$psi_M[i])
  Psi_F2_max = as.numeric(setting_list$psiF2)
  Sigma_A1 = as.numeric(setting_list$sigmaA1)
  Sigma_M1 = as.numeric(setting_list$sigmaM1)
  Sigma_A2 = as.numeric(setting_list$sigmaA2)
  Sigma_M2 = as.numeric(setting_list$sigmaM2)
  e_0 = as.numeric(setting_list$e)
  Q_0 = as.numeric(setting_list$Q)
  
  ## Variables ----
  A1_0 = as.numeric(setting_list$A1_0)
  M1_0 = as.numeric(setting_list$M1_0)
  F1_0 = as.numeric(setting_list$F1_0)
  A2_0 = as.numeric(setting_list$A2_0)
  M2_0 = as.numeric(setting_list$M2_0)
  F2_0 = as.numeric(setting_list$F2_0)
  C_0  = as.numeric(setting_list$C)
  
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
                  gamma_M1.F1_max = Gamma_M1.F1_max,
                  gamma_A2.F2_max = Gamma_A2.F2_max,
                  gamma_M2.F2_max = Gamma_M2.F2_max,
                  gamma_A1.F2_max = Gamma_A1.F2_max,
                  gamma_M1.F2_max = Gamma_M1.F2_max,
                  gamma_A2.F1_max = Gamma_A2.F1_max,
                  gamma_M2.F1_max = Gamma_M2.F1_max,
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
  # Start with initial densities
  results <- data.frame(time = 0, A1 = A1_0, M1 = M1_0, F1 = F1_0,
                        A2 = A2_0, M2 = M2_0, F2 = F2_0, C = C_0,
                        Cycle = NA, Phase = NA)
  
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
        C  = as.numeric(setting_list$C)
        
        state <- c(A1 = A1, 
                   M1 = M1, 
                   F1 = F1, 
                   A2 = A2, 
                   M2 = M2, 
                   F2 = F2, 
                   C  = C)
        
        out <- ode(
          y = state, 
          times = times, 
          func = model, 
          parms = parameters, 
          method = "lsoda", 
          atol = 1e-10,     # Tighten absolute tolerance
          rtol = 1e-8,      # Tighten relative tolerance
          rootfun = rootfunc,
          events = list(func = eventfunc, root = TRUE) # Tells it to use the clamping event
        )
        out_df <- as.data.frame(out) %>%
          mutate(Phase = "growth")

        # Save start and end times of the phase
        df_start <- head(out_df, 1) %>%
          mutate(Type = "start")
        df_end <- tail(out_df, 1) %>%
          mutate(Type = "end")
        
        
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
        
        if (last_phase == "growth" | last_phase == "transconjugant_selection" & Selection_type == "liquid" & time_start != 0){
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
        C  = as.numeric(setting_list$C)
        
        state <- c(A1 = A1, 
                   M1 = M1, 
                   F1 = F1, 
                   A2 = A2, 
                   M2 = M2, 
                   F2 = F2, 
                   C  = C)
        
        out <- ode(
          y = state, 
          times = times, 
          func = model, 
          parms = parameters, 
          method = "lsoda", 
          atol = 1e-10,     # Tighten absolute tolerance
          rtol = 1e-8,      # Tighten relative tolerance
          rootfun = rootfunc,
          events = list(func = eventfunc, root = TRUE) # Tells it to use the clamping event
        )
        out_df <- as.data.frame(out) %>%
          mutate(Phase = "conj")
        
        # Save start and end times of the phase
        df_start <- head(out_df, 1) %>%
          mutate(Type = "start")
        df_end <- tail(out_df, 1) %>%
          mutate(Type = "end")
        
        
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
          C  = as.numeric(setting_list$C)
        } else {
          # For odd cycles, select for F2 transconjugants
          A1 = 0
          M1 = 0
          F1 = 0
          A2 = tail(cycle_df$A2, 1)
          M2 = tail(cycle_df$M2, 1)
          F2 = 0
          C  = as.numeric(setting_list$C)
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
          #total_density = A1 + M1 + A2 + M2
          
          t_out <- tail(cycle_df, 1) %>%
            select(-Phase)
          time_start = tail(t_out$time, 1)
          time_end = time_start + Hours_t_selection
          times <- seq(time_start, time_end, by = Timestep)
          
          # Euler simulation
          state <- c(A1 = A1, 
                     M1 = M1, 
                     F1 = F1, 
                     A2 = A2, 
                     M2 = M2, 
                     F2 = F2, 
                     C  = C)
          
          out <- ode(
            y = state, 
            times = times, 
            func = model, 
            parms = parameters, 
            method = "lsoda", 
            atol = 1e-10,     # Tighten absolute tolerance
            rtol = 1e-8,      # Tighten relative tolerance
            rootfun = rootfunc,
            events = list(func = eventfunc, root = TRUE) # Tells it to use the clamping event
          )
          out_df <- as.data.frame(out)
          
          t_out <- rbind(t_out, out_df)
          
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
          
          
          # Add output to cycle dataframe
          cycle_df <- rbind(cycle_df, t_out)
        } else if (Selection_type == "solid"){
          #### Solid selection: Plate a certain volume of the culture on selective media, then estimate colony size from each cell ----
          time_start = tail(cycle_df, 1)$time
          
          if(Colony_selection == "proportion"){
            ##### Take a proportion of a certain colony number ----
            # Calculate the proportion of each transconjugant type
            A1_por = A1/(A1 + M1)
            A2_por = A2/(A2 + M2)
            M1_por = M1/(A1 + M1)
            M2_por = M2/(A2 + M2)
            
            # Multiple the proportions by x colonies plated
            A1_col = ifelse(is.nan(A1_por), 0, A1_por*Colony_bottleneck)
            A2_col = ifelse(is.nan(A2_por), 0, A2_por*Colony_bottleneck)
            M1_col = ifelse(is.nan(M1_por), 0, M1_por*Colony_bottleneck)
            M2_col = ifelse(is.nan(M2_por), 0, M2_por*Colony_bottleneck)
            
            # Round down to a whole number
            A1_start = ifelse(A1_col < 1, 0, floor(A1_col))/Volume
            A2_start = ifelse(A2_col < 1, 0, floor(A2_col))/Volume
            M1_start = ifelse(M1_col < 1, 0, floor(M1_col))/Volume
            M2_start = ifelse(M2_col < 1, 0, floor(M2_col))/Volume
          } else if(Colony_selection == "volume"){
            ##### Plate out a certain volume ----
            A1_start = ifelse(A1 * T_volume < 1, 0, floor(A1 * T_volume))
            M1_start = ifelse(M1 * T_volume < 1, 0, floor(M1 * T_volume))
            A2_start = ifelse(A2 * T_volume < 1, 0, floor(A2 * T_volume))
            M2_start = ifelse(M2 * T_volume < 1, 0, floor(M2 * T_volume))
          }
          
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
          # Do not dilute strains or add migrants if this is the first phase of the protocol
          
        } else if (last_phase == "transconjugant_selection"){
          # If transitioning from a transconjugant selection phase, do not dilute strains
          # Only relevant for Turner protocol
          
          # Add migrants
          A1 = A1 + A1_migrants
          M1 = M1 + M1_migrants
          F1 = F1 + F1_migrants
          A2 = A2 + A2_migrants
          M2 = M2 + M2_migrants
          F2 = F2 + F2_migrants
          
        } else {
          # Else, dilute strains first
          # Only relevant for Dimitriu protocol
          
          # Add migrants
          A1 = A1*(1 - Immigration_ratio) + A1_migrants*(Immigration_ratio)
          M1 = M1*(1 - Immigration_ratio) + M1_migrants*(Immigration_ratio)
          F1 = F1*(1 - Immigration_ratio) + F1_migrants*(Immigration_ratio)
          A2 = A2*(1 - Immigration_ratio) + A2_migrants*(Immigration_ratio)
          M2 = M2*(1 - Immigration_ratio) + M2_migrants*(Immigration_ratio)
          F2 = F2*(1 - Immigration_ratio) + F2_migrants*(Immigration_ratio)
          
          # Dilute strains
          A1 = ifelse((A1 * Dilution_growth) < Dilution_cutoff, 0, A1 * Dilution_growth)
          M1 = ifelse((M1 * Dilution_growth) < Dilution_cutoff, 0, M1 * Dilution_growth)
          F1 = ifelse((F1 * Dilution_growth) < Dilution_cutoff, 0, F1 * Dilution_growth)
          A2 = ifelse((A2 * Dilution_growth) < Dilution_cutoff, 0, A2 * Dilution_growth)
          M2 = ifelse((M2 * Dilution_growth) < Dilution_cutoff, 0, M2 * Dilution_growth)
          F2 = ifelse((F2 * Dilution_growth) < Dilution_cutoff, 0, F2 * Dilution_growth)
        }
        
        
        # Euler simulation
        C  = as.numeric(setting_list$C)
        
        state <- c(A1 = A1, 
                   M1 = M1, 
                   F1 = F1, 
                   A2 = A2, 
                   M2 = M2, 
                   F2 = F2, 
                   C  = C)
        
        out <- ode(
          y = state, 
          times = times, 
          func = model, 
          parms = parameters, 
          method = "lsoda", 
          atol = 1e-10,     # Tighten absolute tolerance
          rtol = 1e-8,      # Tighten relative tolerance
          rootfun = rootfunc,
          events = list(func = eventfunc, root = TRUE) # Tells it to use the clamping event
        )
        out_df <- as.data.frame(out) %>%
          mutate(Phase = "immigration")
        
        # Save start and end times of the phase
        df_start <- head(out_df, 1) %>%
          mutate(Type = "start")
        df_end <- tail(out_df, 1) %>%
          mutate(Type = "end")
        
        
        # Add output to cycle dataframe
        cycle_df <- rbind(cycle_df, out_df)
      } 
      # End simulation if ancestor or mutant density is 0
      final_A = tail(cycle_df$A1, 1) + tail(cycle_df$A2, 1)
      final_M = tail(cycle_df$M1, 1) + tail(cycle_df$M2, 1)
      
      if(final_A == 0 | final_M == 0){
        break
      }
    }
    ### Cycle wrap up ----
    cycle_df <- cycle_df %>%
      mutate(Cycle = c) %>%
      slice(-1)
    
    # Save results to main dataframes
    results <- rbind(results, cycle_df) %>%
      drop_na()
    
    # End simulation if ancestor or mutant density is 0
    final_A = tail(results$A1, 1) + tail(results$A2, 1)
    final_M = tail(results$M1, 1) + tail(results$M2, 1)
    
    if(final_A == 0 | final_M == 0){
      break
    }
    
  }
  # Save invasion information
  sw <- sweep_param[i,] %>%
    mutate(end_time = tail(results,1)$time,
           A1 = tail(results,1)$A1,
           M1 = tail(results,1)$M1,
           A2 = tail(results,1)$A2,
           M2 = tail(results,1)$M2)
  sw
  #sweep_out <- rbind(sweep_out,sw)
  
        }

stopCluster(cluster)

# Save output file
write_csv(sweep_out, paste0(output_folder, "/" , ps, "_out.csv"))
