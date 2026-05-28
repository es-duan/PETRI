# Process plating growth curve data 

# Load in packages ----
library(tidyverse)
library(broom)
library(argparse)
library(jsonlite)

# Set arguments parser inputs ----
parser <- ArgumentParser()
parser$add_argument("-c","--colors", help = "JSON string of plot colors")

# Parse arguments
args <- parser$parse_args()

# Load global variables ----
## Colors ----
plot_colors <- fromJSON(args$colors)
p_Anc <- plot_colors[["p_Anc"]]
p_Mut <- plot_colors[["p_Mut"]]
p_F <- plot_colors[["p_F"]]

## Retrieve ggplot theme ----
source("src/ggplot_theme.R")

# Read in plating file ----
plating <- read_csv("input_data/experimental_data/2026-04-02_growth_rate_plating.csv")
output_dir <- "results/phenotyping/growth_rate/plating"

# Select relevant data and recalculate density ----
plating2 <- plating %>%
  filter(Strain %in% c("TR44", "TR84", "TR85")) %>%
  mutate(Density = Count * (1000/Volume_plated) * 10^Plate_Dilution) %>%
  mutate(Strain = case_when(Strain == "TR44" ~ "F",
                            Strain == "TR84" ~ "Anc",
                            Strain == "TR85" ~ "Mut"))

# Plot density over time ----
p1 <- ggplot(data = plating2,
       mapping = aes(Time, Density, color = Strain)) +
  stat_summary(fun.data="mean_se",geom="errorbar",width=0.2) +
  stat_summary(fun = "mean", geom = "line") +
  stat_summary(fun = "mean", geom = "point") +
  scale_y_continuous(trans = "log10") +
  scale_color_manual(values = c("F" = p_F,
                                "Anc" = p_Anc,
                                "Mut" = p_Mut)) +
  fig_aes
ggsave(paste(output_dir, "growth_density_av.pdf", sep = "/"),
       p1, width = 6, height = 4, units = "in")


# Calculate growth rate over 4 hour window ----
growth <- plating2 %>%
  select(Strain, Replicate, Time, Density) %>%
  pivot_wider(names_from = Time, values_from = Density, names_prefix = "Dens_") %>%
  mutate(growth_rate = log(Dens_4/Dens_0)/4)

## Growth rate plot ----
p2 <- ggplot(data = growth,
             mapping = aes(Strain, growth_rate, color = Strain)) +
  stat_summary(fun.data="mean_se",geom="errorbar",width=0.2) +
  stat_summary(fun = "mean", geom = "line") +
  stat_summary(fun = "mean", geom = "point") +
  scale_color_manual(values = c("F" = p_F,
                                "Anc" = p_Anc,
                                "Mut" = p_Mut)) +
  fig_aes
ggsave(paste(output_dir, "growth_rate_av.pdf", sep = "/"),
       p2, width = 5, height = 4, units = "in")

# Average data ----
gr_out <- growth %>%
  group_by(Strain) %>%
  summarise(gr_mean = mean(growth_rate),
            gr_sd = sd(growth_rate),
            n = n()) %>%
  ungroup() %>%
  mutate(gr_se = gr_sd/sqrt(n))
write_csv(gr_out, paste(output_dir, "growthrate0-4_av.csv", sep = "/"))

## Statistics ----
gr_tt <- t.test(filter(growth, Strain == "Anc")$growth_rate, filter(growth, Strain == "Mut")$growth_rate, 
                alternative = "two.sided",
                paired = TRUE,
                var.equal = FALSE,
                conf.level = 0.95)
write_csv(tidy(gr_tt), paste(output_dir, "growthrate_tt.csv", sep = "/"))

