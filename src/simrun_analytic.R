# Script to run the analytic model

# Load packages ----
library(deSolve)
library(tidyverse)
library(argparse)

# Set arguments parser inputs ----
parser <- ArgumentParser()
parser$add_argument("-t","--treatment", type = "character", help = "Specify Treatment ID")
args <- parser$parse_args()

treatment <- args$treatment

# Set parameters ----

## Load parameter file ----
treatment_csv <- read.csv("input_data/analytic_parameters.csv", header = F)

## Transpose data to make it easier for reading in parameters ----
tr_colnames <- treatment_csv[[1]]
treatment_csv <- as.data.frame(t(treatment_csv[,-1]))
colnames(treatment_csv) <- tr_colnames

## Parameters
Gamma_A.F = 0.01
Gamma_M.F_change = 0.6
Psi_A = 1
Psi_M_cost = 0.1
Psi_F = 1.1
Sigma = 0.25
Delta = 0.1
Alpha = 0

## Variables ----
A_0 = 1e9
M_0 = 1
F_0 = 1e8

# Express densities as proportions
a_0 <- A_0/(A_0 + M_0 + F_0)
m_0 <- M_0/(A_0 + M_0 + F_0)
f_0 <- F_0/(A_0 + M_0 + F_0)

# Define the model ----
model <- function(time, state, parameters){
  with(as.list(c(state, parameters)), {
    # Set model state
    a <- state["a"]
    m <- state["m"]
    f <- state["f"]
    
    # Differential equations
    da <- (psi_A*(1-sigma)*a) + (gamma_A.F*a*f) - (delta*a) - (a*((psi_A-delta)*a + (psi_A - psi_M_cost - delta)*m - alpha*f))
    dm <- ((psi_A - psi_M_cost)*(1-sigma)*m) + ((gamma_A.F + gamma_M.F_change)*m*f) - (delta*m) - (m*((psi_A-delta)*a + (psi_A - psi_M_cost - delta)*m - alpha*f))
    df <- -(delta*f) + (psi_A*sigma*a) + ((psi_A - psi_M_cost)*sigma*m) - (gamma_A.F*a*f) - ((gamma_A.F + gamma_M.F_change)*m*f) - (f*((psi_A-delta)*a + (psi_A - psi_M_cost - delta)*m - alpha*f))
    
    list(c(da, dm, df))
  })
}

# Test run
state <- c(a = a_0,
           m = m_0,
           f = f_0)

parameters <- c(gamma_A.F = Gamma_A.F,
                gamma_M.F_change = Gamma_M.F_change,
                psi_A = Psi_A,
                psi_M_cost = Psi_M_cost,
                psi_F = Psi_F,
                sigma = Sigma,
                delta = Delta,
                alpha = Alpha)

times <- seq(0, 600, 0.1)
out <- ode(y = state, times = times, func = model, parms = parameters, method = "euler")
out_df <- as.data.frame(out)

# Plot data to visualize
data_long <- out_df %>%
  pivot_longer(cols = c(a, m, f),
               names_to = "Cell_type",
               values_to = "Proportion")

ggplot(data_long, aes(time, Proportion, color = Cell_type)) +
  geom_line()
