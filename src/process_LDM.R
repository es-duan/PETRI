# Process LDM data

# Load packages ----
library(tidyverse)
library(openxlsx)
library(broom)
library(argparse)
library(jsonlite)

# Set arguments parser inputs ----
parser <- ArgumentParser()
parser$add_argument("-c","--colors", help = "JSON string of plot colors")
parser$add_argument("-o","--points", help = "JSON string of point aesthetics")

# Parse arguments
args <- parser$parse_args()

# Load global variables ----
## Colors ----
plot_colors <- fromJSON(args$colors)
p_Anc <- plot_colors[["p_Anc"]]
p_Mut <- plot_colors[["p_Mut"]]
p_F <- plot_colors[["p_F"]]

## Points ----
plot_points <- jsonlite::fromJSON(args$points)
exp_point_size <- plot_points[["exp_point_size"]]
sh_Anc <- plot_points[["sh_Anc"]]
sh_Mut <- plot_points[["sh_Mut"]]
sh_F <- plot_points[["sh_F"]]

## Retrieve ggplot theme ----
source("src/ggplot_theme.R")

# Read in data ----
density <- read.xlsx("input_data/experimental_data/2026-02-18_compiled_LDM_plating.xlsx",
                     sheet = "2026-02-18_compiled_LDM_plating", detectDates = TRUE)
turbidity <- read.xlsx("input_data/experimental_data/2026-02-18_compiled_LDM_plating.xlsx",
                       sheet = "turbidity", detectDates = TRUE)

# Calculate density averages ----
density2 <- density %>%
  mutate(Density = Counts * (1000/Volume.plated) * 10^Dilution) %>%
  group_by(Day, Strain, Plate.cell.type, Time) %>%
  summarize(Density_mean = mean(Density),
            Density_sd = sd(Density),
            n = n()) %>%
  ungroup() %>%
  mutate(Density_se = Density_sd/sqrt(n))

## Plot densities ----
density_plot <- density2 %>%
  mutate(Strain = case_when(Strain == "84" & Plate.cell.type == "Donor" ~ "Anc",
                              Strain == "85" & Plate.cell.type == "Donor" ~ "Mut",
                              Plate.cell.type == "Recipient" ~ "F"))

p1 <- ggplot(data = density_plot,
       mapping = aes(Time, Density_mean, color = Strain, shape = Strain)) +
  geom_line() +
  geom_errorbar(data = density_plot,
                mapping = aes(ymin = Density_mean - Density_se, ymax = Density_mean + Density_se),
                width = 0.1) +
  geom_point(size = exp_point_size) +
  scale_color_manual(values = c("F" = p_F,
                                "Anc" = p_Anc,
                                "Mut" = p_Mut)) +
  scale_shape_manual(values = c("F" = sh_F,
                                "Anc" = sh_Anc,
                                "Mut" = sh_Mut)) +
  facet_grid(~Day) +
  fig_aes

ggsave("results/phenotyping/LDM_conjugation/LDM_density.pdf",
       p1, width = 10, height = 4, units = "in")

# Merge with turbidity and calculate LDM ----
LDM <- density2 %>%
  select(-Density_sd, -Density_se, -n) %>%
  mutate(Time = case_when(Time == 0 ~ "0",
                          TRUE ~"final")) %>%
  pivot_wider(names_from = c("Plate.cell.type", "Time"), values_from = Density_mean) %>%
  left_join(select(turbidity, Day, Strain, Time_h, p0), by = c("Day", "Strain")) %>%
  mutate(f = 10) %>% # 100 uL matings
  mutate(LDM = ((f)*(1/Time_h)*(-log(p0))*((log(Donor_final*Recipient_final)-log(Donor_0*Recipient_0))/((Donor_final*Recipient_final)-(Donor_0*Recipient_0))))) %>%
  mutate(Strain = case_when(Strain == 84 ~ "Anc",
                            Strain == 85 ~ "Mut"))

LDM_av <- LDM %>%
  group_by(Strain) %>%
  summarize(Conjugation_rate_mean = mean(LDM),
            Conjugation_rate_sd = sd(LDM),
            n = n()) %>%
  ungroup() %>%
  mutate(Conjugation_rate_se = Conjugation_rate_sd/sqrt(n))

## Plot final conjugation rates ----
p2 <- ggplot(data = LDM_av,
       mapping = aes(Strain, Conjugation_rate_mean, color = Strain, shape = Strain)) +
  geom_line() +
  geom_errorbar(data = LDM_av,
                mapping = aes(ymin = Conjugation_rate_mean - Conjugation_rate_se, ymax = Conjugation_rate_mean + Conjugation_rate_se),
                width = 0.1) +
  geom_point(size = exp_point_size) +
  scale_color_manual(values = c("F" = p_F,
                                "Anc" = p_Anc,
                                "Mut" = p_Mut)) +
  scale_shape_manual(values = c("F" = sh_F,
                                "Anc" = sh_Anc,
                                "Mut" = sh_Mut)) +
  scale_y_continuous(trans = "log10") +
  fig_aes

ggsave("results/phenotyping/LDM_conjugation/LDM_conj_rates.pdf",
       p2, width = 4, height = 3, units = "in")

# Statistical analysis ----
conj_tt <- t.test(log10(filter(LDM, Strain == "Anc")$LDM),
                  log10(filter(LDM, Strain == "Mut")$LDM), 
                  alternative = "two.sided",
                  var.equal = FALSE,
                  paired = TRUE,
                  conf.level = 0.95)
write_csv(tidy(conj_tt), "results/phenotyping/LDM_conjugation/LDM_conjugation_tt.csv")

# Save final file ----
write_csv(LDM_av, "results/phenotyping/LDM_conjugation/LDM_conjugation_av.csv")

