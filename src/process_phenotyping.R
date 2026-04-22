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

## Read in growth rate data ----
# Update script to read file once growth is finalized
OD_growth <- data.frame("Genotype" = c("Anc", "Mut", "F"),
                        "Growth_rate_mean" = c(0.476, 0.623, 0.607),
                        "Growth_rate_sd" = c(0.0323, 0.0440, 0.0258),
                        "n" = rep("n", 3),
                        "Growth_rate_se" = c(0.0186, 0.0254, 0.0149))

## Conjugation rate dataframe ----
LDM <- data.frame("Genotype" = c(rep("Anc",3), rep("Mut", 3)),
                  "Conjugation_rate" = c(1.56e-12, 1.16e-12, 4.60e-12,
                                         3.53e-14, 9.12e-15, 5.38e-15),
                  "Day" = as.Date(rep(c("2025-07-23","2025-08-25","2025-09-11"), 2)))

# Summarize LDM data ----
LDM_av <- LDM %>%
  group_by(Genotype) %>%
  summarise(Conjugation_rate_mean = mean(Conjugation_rate),
            Conjugation_rate_sd = sd(Conjugation_rate),
            n = n()) %>%
  ungroup() %>%
  mutate(Conjugation_rate_se = Conjugation_rate_sd/sqrt(n))

# Merge growth and conjugation data ----
phenotyping <- left_join(OD_growth, LDM_av, by = "Genotype")

# Statistics ----
## Growth rate ----
growth_tt <- t.test(filter(OD_growth, Genotype == "Anc")$Growth_rate,
                    filter(OD_growth, Genotype == "Mut")$Growth_rate, 
                 alternative = "two.sided",
                 var.equal = FALSE,
                 paired = TRUE,
                 conf.level = 0.95)

## Conjugation rate ----
conj_tt <- t.test(log10(filter(LDM, Genotype == "Anc")$Conjugation_rate),
                  log10(filter(LDM, Genotype == "Mut")$Conjugation_rate), 
                  alternative = "two.sided",
                  var.equal = FALSE,
                  paired = TRUE,
                  conf.level = 0.95)


# Export data ----
write_csv(LDM_phenotype, "results/phenotyping/LDM_phenotyping.csv")
write_csv(phenotyping, "results/phenotyping/phenotyping_av.csv")

write_csv(tidy(growth_tt), "results/phenotyping/LDM_growth_tt.csv")
write_csv(tidy(conj_tt), "results/phenotyping/LDM_conjugation_tt.csv")

