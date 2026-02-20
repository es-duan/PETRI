# Process and plot strain phenotyping data

# Load packages ----
library(tidyverse)
library(broom)

# Read in data ----
LDM_plating <- read_csv("input_data/experimental_data/2026-02-18_compiled_LDM_plating.csv")[1:10]

# Process data ----
## Change names for plotting ----
LDM_plating2 <- LDM_plating 

## Average data ----
LDM_plating_av <- LDM_plating2 %>%
  group_by(Day, Strain, Time, Plate.cell.type) %>%
  summarize(Density_mean = mean(Density),
            Density_sd = sd(Density),
            n = n()) %>%
  ungroup() %>%
  mutate(Density_se = Density_sd/sqrt(n))

## Calculate growth rates ----
LDM_growth <- LDM_plating_av %>%
  select(Day, Strain, Time, Plate.cell.type, Density_mean) %>%
  pivot_wider(names_from = Time, values_from = Density_mean, names_prefix = "Density_") %>%
  mutate(Time = ifelse(!is.na(Density_1), 1, 2),
         Density_t = ifelse(!is.na(Density_1), Density_1, Density_2)) %>%
  select(-Density_1, -Density_2) %>%
  mutate(Growth_rate = log(Density_t/Density_0)/Time)

LDM_growth_av <- LDM_growth %>%
  group_by(Strain, Plate.cell.type) %>%
  summarize(Growth_rate_mean = mean(Growth_rate),
            Growth_rate_sd = sd(Growth_rate),
            n = n()) %>%
  ungroup() %>%
  mutate(Growth_rate_se = Growth_rate_sd/sqrt(n))

## Conjugation rate dataframe ----
LDM <- data.frame("Genotype" = c(rep("Anc",3), rep("Mut", 3)),
                  "Conjugation_rate" = c(1.56e-12, 1.16e-12, 4.60e-12,
                                         3.53e-14, 9.12e-15, 5.38e-15),
                  "Day" = as.Date(rep(c("2025-07-23","2025-08-25","2025-09-11"), 2)))

# Combine data ----
LDM_phenotype <- LDM_growth %>%
  mutate(Genotype = case_when(Strain == 84 & Plate.cell.type == "Donor" ~ "Anc",
                              Strain == 85 & Plate.cell.type == "Donor" ~ "Mut",
                              Plate.cell.type == "Recipient" ~ "F")) %>%
  filter(Genotype != "F") %>%
  select(Day, Genotype, Growth_rate) %>%
  left_join(LDM, by = c("Day", "Genotype"))

LDM_phenotype_av <- LDM_phenotype %>%
  group_by(Genotype) %>%
  summarise(Growth_rate_mean = mean(Growth_rate),
            Growth_rate_sd = sd(Growth_rate),
            Conjugation_rate_mean = mean(Conjugation_rate),
            Conjugation_rate_sd = sd(Conjugation_rate),
            n = n()) %>%
  ungroup() %>%
  mutate(Growth_rate_se = Growth_rate_sd/sqrt(n),
         Conjugation_rate_se = Conjugation_rate_sd/sqrt(n))

# Statistics ----
## Growth rate ----
growth_tt <- t.test(filter(LDM_phenotype, Genotype == "Anc")$Growth_rate,
                    filter(LDM_phenotype, Genotype == "Mut")$Growth_rate, 
                 alternative = "two.sided",
                 var.equal = FALSE,
                 paired = TRUE,
                 conf.level = 0.95)

## Conjugation rate ----
conj_tt <- t.test(log10(filter(LDM_phenotype, Genotype == "Anc")$Conjugation_rate),
                  log10(filter(LDM_phenotype, Genotype == "Mut")$Conjugation_rate), 
                  alternative = "two.sided",
                  var.equal = FALSE,
                  paired = TRUE,
                  conf.level = 0.95)


# Export data ----
write_csv(LDM_phenotype, "results/phenotyping/LDM_phenotyping.csv")
write_csv(LDM_phenotype_av, "results/phenotyping/LDM_phenotyping_av.csv")

write_csv(tidy(growth_tt), "results/phenotyping/LDM_growth_tt.csv")
write_csv(tidy(conj_tt), "results/phenotyping/LDM_conjugation_tt.csv")

