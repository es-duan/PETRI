# Process experimental data

# Load packages ----
library(tidyverse)
library(broom)

# Data processing/experimental parameters ----
tselect_t = 72 # colony growth time of 72 hr
tselect_d = 1e7 # scale plate scraping estimates to a density of 1e7
limit_est = 0.5 # use estimate of 0.5 for plates below the detection limit

## Relevel factors ----
host_levels <- c("rifR", "nalR")
ph_levels <- c("Anc", "Mut", "F")
phase_levels <- c("growth", "conj", "tselect")

# Read in data ----
plating_HFC <- read_csv("input_data/experimental_data/2026-01-27_HFC_plating.csv")
plating_LFC <- read.csv("input_data/experimental_data/2026-01-27_LFC_plating.csv")

extinction_out_av <- read_csv("results/experimental_validation/extinction_cell_counts/extinction_size_av.csv") %>%
  select(Host, Phenotype, Plate_type, Extinction_mean, Cells_colony_mean)

# Specify output directory
output_dir <- "results/experimental_validation"

# Process HFC data ----
## Clean up data, apply extinctions, add metadata ----
plating_HFC2 <- plating_HFC %>%
  filter(!is.na(Count)) %>%
  mutate(Replicate = as.character(Replicate)) %>%
  mutate(Host = ifelse(str_detect(Plate_type, "Rif"), "rifR", "nalR")) %>% # Add host info
  left_join(extinction_out_av, by = c("Host", "Phenotype", "Plate_type")) %>% # Add extinction data
  mutate(Count_c = Count/(1-Extinction_mean)) %>% # Apply extinction
  mutate(Density_c = Count_c * (1000/Volume_plated) * (10^Plate_Dilution)) # Recalculate density


## Calculate colony size estimates ----
HFC_t <- plating_HFC2 %>%
  filter(Time == 5 & Host == "nalR" & Phenotype %in% c("Anc", "Mut")) %>%
  mutate(Colony_density_est = (Count_c*Cells_colony_mean)/5) %>% # Assume resuspending in 5 mL
  select(Experiment, Treatment, 
         Time, Replicate, Phenotype, Colony_density_est) %>%
  pivot_wider(names_from = Phenotype, values_from = Colony_density_est) %>%
  mutate(Total_density = Anc + Mut,
         Dil_factor = Total_density/tselect_d) %>%
  mutate(Anc_end = Anc/Dil_factor,
         Mut_end = Mut/Dil_factor) %>%
  select(Experiment, Treatment,
         Time, Replicate, Anc_end, Mut_end) %>%
  pivot_longer(c(Anc_end, Mut_end), names_to = "Phenotype", values_to = "Density_c") %>%
  mutate(Time = 5 + tselect_t,
         Phenotype = str_replace(Phenotype, "_end", ""),
         Host = "nalR",
         Phase = "tselect",
         Count_type = "estimate")

## Merge data ----
HFC_out <- plating_HFC2 %>%
  bind_rows(HFC_t) %>%
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
write_csv(HFC_out, paste(output_dir, "HFC", "HFC_plating_processed.csv", sep = "/"))
write_csv(HFC_out_av, paste(output_dir, "HFC", "HFC_plating_processed_av.csv", sep = "/"))
write_csv(HFC_f, paste(output_dir, "HFC", "HFC_frequency_processed.csv", sep = "/"))
write_csv(HFC_f_av, paste(output_dir, "HFC", "HFC_frequency_processed_av.csv", sep = "/"))
write_csv(tidy(HFC_tt), paste(output_dir, "HFC", "HFC_ttest.csv", sep = "/"))

# Process LFC data ----
## Clean up data, apply extinctions, add metadata ----
plating_LFC2 <- plating_LFC %>%
  filter(!is.na(Count)) %>%
  # Select first 3 replicates only
  filter(Replicate %in% c(1,2,3)) %>%
  mutate(Replicate = as.character(Replicate)) %>%
  mutate(Host = ifelse(str_detect(Plate_type, "Rif"), "rifR", "nalR")) %>% # Add host info
  left_join(extinction_out_av, by = c("Host", "Phenotype", "Plate_type")) %>% # Add extinction data
  mutate(Count_c = Count/(1-Extinction_mean)) %>% # Apply extinction
  mutate(Count_c = ifelse(Count_type == "below_limit", limit_est, Count_c)) %>% # Add estimate for points below detection limit
  mutate(Density_c = Count_c * (1000/Volume_plated) * (10^Plate_Dilution)) # Recalculate density


## Calculate colony size estimates ----
LFC_t <- plating_LFC2 %>%
  filter(Time == 173 & Host == "nalR" & Phenotype %in% c("Anc", "Mut")) %>%
  mutate(Colony_density_est = (Count_c*Cells_colony_mean)/5) %>% # Assume resuspending in 5 mL
  select(Experiment, Treatment, 
         Time, Replicate, Phenotype, Colony_density_est) %>%
  pivot_wider(names_from = Phenotype, values_from = Colony_density_est) %>%
  mutate(Total_density = Anc + Mut,
         Dil_factor = Total_density/tselect_d) %>%
  mutate(Anc_end = Anc/Dil_factor,
         Mut_end = Mut/Dil_factor) %>%
  select(Experiment, Treatment,
         Time, Replicate, Anc_end, Mut_end) %>%
  pivot_longer(c(Anc_end, Mut_end), names_to = "Phenotype", values_to = "Density_c") %>%
  mutate(Time = 173 + tselect_t,
         Phenotype = str_replace(Phenotype, "_end", ""),
         Host = "nalR",
         Phase = "tselect",
         Count_type = "estimate")

## Estimate T168 conj dilution ----
LFC_c <- plating_LFC2 %>%
  filter(Time == 168 & Phase == "growth") %>%
  mutate(Density_c = Density_c/100) %>% # 100-fold dilution into conjugation phase
  mutate(Phase = "conj",
         Count_type = "dilution_est")

## Merge data ----
LFC_out <- plating_LFC2 %>%
  bind_rows(LFC_t, LFC_c) %>%
  select(Experiment, Treatment, Phase, Time, Replicate, Phenotype, Host,
         Count, Count_c, Density_c, Count_type) %>%
  mutate(Phase = factor(Phase, levels = phase_levels),
         Host = factor(Host, levels = host_levels),
         Phenotype = factor(Phenotype, levels = ph_levels))

## Calculate frequencies ----
LFC_f <- LFC_out %>%
  filter(Phenotype %in% c("Anc", "Mut")) %>%
  filter(!(Time == 173 & Host == "rifR")) %>%
  filter(!(Time == 168 & Phase == "conj")) %>%
  select(-Count, -Count_c, -Count_type) %>%
  pivot_wider(names_from = Phenotype, values_from = Density_c) %>%
  mutate(Density_p = Anc + Mut,
         Anc_freq = Anc/Density_p,
         Mut_freq = Mut/Density_p) %>%
  pivot_longer(c(Anc_freq, Mut_freq), names_to = "Phenotype", values_to = "Frequency") %>%
  mutate(Phenotype = str_replace(Phenotype, "_freq", "")) %>%
  select(-Anc, -Mut)

## Determine count type for averages ----
LFC_ct_times <- LFC_out %>% 
  filter(Count_type != "count") %>%
  select(Phase, Time, Phenotype, Host, Count_type) %>%
  distinct()

LFC_ct <- LFC_out %>%
  select(Phase, Time, Phenotype, Host) %>%
  distinct() %>%
  left_join(LFC_ct_times, by = c("Time", "Phase", "Phenotype", "Host")) %>%
  mutate(Count_type = ifelse(is.na(Count_type), "count", Count_type))


## Average data ----
LFC_out_av <- LFC_out %>%
  group_by(Experiment, Treatment, Phase, Time, Phenotype, Host) %>%
  summarize(Density_mean = mean(Density_c),
            Density_sd = sd(Density_c),
            n = n()) %>%
  ungroup() %>%
  mutate(Density_se = Density_sd/sqrt(n)) %>%
  left_join(LFC_ct, by = c("Time", "Phase", "Phenotype", "Host"))

LFC_f_av <- LFC_f %>%
  group_by(Experiment, Treatment, Phase, Time, Phenotype, Host) %>%
  summarize(Frequency_mean = mean(Frequency),
            Frequency_sd = sd(Frequency),
            n = n()) %>%
  ungroup() %>%
  mutate(Frequency_se = Frequency_sd/sqrt(n)) %>%
  left_join(LFC_ct, by = c("Time", "Phase", "Phenotype", "Host"))

## Statistical analysis ----
LFC_stat <- LFC_f %>%
  filter(Time %in% c(0,173+tselect_t)) %>%
  filter(Phenotype == "Mut") %>%
  select(Time, Phenotype, Frequency, Replicate) %>%
  pivot_wider(names_from = Time, values_from = Frequency,
              names_prefix = "Freq_T") %>%
  mutate(Freq_T0a = asin(sqrt(Freq_T0)),
         !!paste0("Freq_T",173+tselect_t,"a") := asin(sqrt(!!sym(paste0("Freq_T",173+tselect_t))))) # arcsin transformation

LFC_tt <- t.test(LFC_stat[[paste0("Freq_T",173+tselect_t,"a")]], LFC_stat$Freq_T0a, 
       alternative = "greater",
       paired = TRUE,
       var.equal = FALSE,
       conf.level = 0.95)

## Save files ----
write_csv(LFC_out, paste(output_dir, "LFC", "LFC_plating_processed.csv", sep = "/"))
write_csv(LFC_out_av, paste(output_dir, "LFC", "LFC_plating_processed_av.csv", sep = "/"))
write_csv(LFC_f, paste(output_dir, "LFC", "LFC_frequency_processed.csv", sep = "/"))
write_csv(LFC_f_av, paste(output_dir, "LFC", "LFC_frequency_processed_av.csv", sep = "/"))
write_csv(tidy(LFC_tt), paste(output_dir, "LFC", "LFC_ttest.csv", sep = "/"))

