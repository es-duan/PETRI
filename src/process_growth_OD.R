# Process OD growth curve data 

# Load in packages ----
library(tidyverse)
library(openxlsx)
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

# Read in OD file ----
file_name <- "input_data/experimental_data/2026-04-16_OD_growth_curve.xlsx"
start_row <- 53 #Row of Cycle data
end_row <- 109 #Last well
loop_time <- 5 #Time interval of each kinetic loop
output_dir <- "results/phenotyping/growth_rate/OD"

growth <- read.xlsx(file_name,sheet = "Result sheet",
                    rows = start_row:end_row)[-1,][-1,]

format <- read.xlsx(file_name,sheet = "Metadata")

# Reformat to long form ----
growth2 <- growth %>%
  rename(Well = Cycle.Nr.) %>%
  pivot_longer(!Well, names_to = "Cycle", values_to = "OD600") %>%
  select(Well,OD600,Cycle) %>%
  left_join(format, by = "Well") %>%
  filter(Strain %in% c("TR44", "TR84", "TR85")) %>%
  mutate(Strain = case_when(Strain == "TR44" ~ "F",
                            Strain == "TR84" ~ "Anc",
                            Strain == "TR85" ~ "Mut")) %>%
  mutate_at(c("Replicate", "Rep2"), as.character) %>%
  mutate(Time = round(((as.numeric(Cycle)*loop_time))/60,1)) %>%
  # Look at growth curves for first 20 hours
  filter(Time <= 20) %>%
  # Second dilution replicate produced the most consistent results
  filter(Rep2 == "2") %>%
  select(-Rep2)

## Plot all OD data ----
### Smooth plot ----
p1 <- ggplot(data = growth2,
       mapping = aes(Time, OD600, color = Strain, linetype = Replicate)) +
  geom_line() +
  facet_grid(~Strain) +
  scale_color_manual(values = c("F" = p_F,
                                "Anc" = p_Anc,
                                "Mut" = p_Mut)) +
  fig_aes
ggsave(paste(output_dir, "growthcurve_allrep.pdf", sep = "/"),
       p1, width = 10, height = 4, units = "in")

### Summarized plot with points at interval used for max gr calculation ----
p2 <- ggplot(data = filter(growth2, Time %in% c(0.1, seq(1, 20, 0.5))),
       mapping = aes(Time, OD600, color = Strain)) +
  stat_summary(fun.data="mean_se",geom="errorbar",width=0.2) +
  stat_summary(fun = "mean", geom = "line") +
  stat_summary(fun = "mean", geom = "point") +
  #scale_y_continuous(trans = "log10") +
  scale_color_manual(values = c("F" = p_F,
                                "Anc" = p_Anc,
                                "Mut" = p_Mut)) +
  fig_aes
ggsave(paste(output_dir, "growthcurve_av.pdf", sep = "/"),
       p2, width = 8, height = 4, units = "in")
saveRDS(p2, paste(output_dir, "growthcurve_av.rds", sep = "/"))

# Calculate growth rates over 30 min intervals ----
growth_rate <- growth2 %>%
  select(Strain, Replicate, Time, OD600) %>%
  arrange(Strain, Replicate, Time) %>%
  filter(Time %in% c(0.1, seq(1, 20, 0.5))) %>%
  group_by(Strain, Replicate) %>%
  mutate(growth_rate = log(OD600/lag(OD600))/(Time - lag(Time))) %>%
  ungroup()

### Plot growth rates ----
p3 <- ggplot(growth_rate,
       aes(Time, growth_rate, color = Strain)) +
  stat_summary(fun.data="mean_se",geom="errorbar",width=0.2, alpha = 0.75) +
  stat_summary(fun = "mean", geom = "line") +
  stat_summary(fun = "mean", geom = "point") +
  ylab("Growth Rate") +
  #scale_y_continuous(trans = "log10") +
  scale_color_manual(values = c("F" = p_F,
                                "Anc" = p_Anc,
                                "Mut" = p_Mut)) +
  fig_aes
ggsave(paste(output_dir, "growthrate_av.pdf", sep = "/"),
       p3, width = 8, height = 4, units = "in")
saveRDS(p3, paste(output_dir, "growthrate_av.rds", sep = "/"))

# Select max growth rate ----
gr_max <- growth_rate %>%
  group_by(Strain, Replicate) %>%
  slice_max(growth_rate) %>%
  ungroup() %>%
  rename(growth_max = growth_rate)

### Plot max growth rates ----
p4 <- ggplot(gr_max, aes(Strain, growth_max, color = Strain)) +
  stat_summary(fun.data="mean_se",geom="errorbar",
               width = 0.1, alpha = 0.75) +
  stat_summary(fun = "mean", geom = "point") +
  scale_y_continuous(limits = c(0.4, 0.7),
                     name = "max growth rate") +
  scale_color_manual(values = c("F" = p_F,
                                "Anc" = p_Anc,
                                "Mut" = p_Mut)) +
  fig_aes
ggsave(paste(output_dir, "growthratemax_av.pdf", sep = "/"),
       p4, width = 8, height = 4, units = "in")
saveRDS(p4, paste(output_dir, "growthratemax_av.rds", sep = "/"))

# Summarize data, output data and statistics ----
gr_out <- gr_max %>%
  group_by(Strain) %>%
  summarise(max_gr_mean = mean(growth_max),
            max_gr_sd = sd(growth_max),
            n = n()) %>%
  ungroup() %>%
  mutate(max_gr_se = max_gr_sd/sqrt(n))
write_csv(gr_out, paste(output_dir, "growthratemax_av.csv", sep = "/"))

## Statistics ----
gr_tt <- t.test(filter(gr_max, Strain == "Anc")$growth_max, filter(gr_max, Strain == "Mut")$growth_max, 
                   alternative = "two.sided",
                   paired = TRUE,
                   var.equal = FALSE,
                   conf.level = 0.95)
write_csv(tidy(gr_tt), paste(output_dir, "growthratemax_tt.csv", sep = "/"))
