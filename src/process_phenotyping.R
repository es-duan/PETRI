# Process and plot strain phenotyping data

# Load packages ----
library(tidyverse)
library(broom)

# Read in data ----
gr_out <- read_csv("results/phenotyping/growth_rate/growthratemax_av.csv")

# Rename data for merging ----
OD_growth <- gr_out %>%
  rename(Genotype = Strain,
         Growth_rate_mean = max_gr_mean,
         Growth_rate_sd = max_gr_sd,
         Growth_rate_se = max_gr_se)

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
## Conjugation rate ----
conj_tt <- t.test(log10(filter(LDM, Genotype == "Anc")$Conjugation_rate),
                  log10(filter(LDM, Genotype == "Mut")$Conjugation_rate), 
                  alternative = "two.sided",
                  var.equal = FALSE,
                  paired = TRUE,
                  conf.level = 0.95)


# Export data ----
write_csv(phenotyping, "results/phenotyping/phenotyping_av.csv")

write_csv(tidy(conj_tt), "results/phenotyping/LDM_conjugation_tt.csv")

