# Generate Fig psweep

# Load packages ----
library(tidyverse)
library(patchwork)
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
plot_colors <- jsonlite::fromJSON(args$colors)
p_Anc <- plot_colors[["p_Anc"]]
p_Mut <- plot_colors[["p_Mut"]]

## Points ----
plot_points <- jsonlite::fromJSON(args$points)
psweep_point_size <- plot_points[["psweep_point_size"]]
sh_Anc <- plot_points[["sh_Anc"]]
sh_Mut <- plot_points[["sh_Mut"]]

## Common aesthetics ----
anc_name <- "X"
mut_name <- "Y"
label_size <- 5
point_color <- "black" # overriding colors

# Read in figures ----
anc_HFC <- readRDS("results/parameter_sweeps/pHFC_S.pB10/pHFC_S.pB10_inv_change_plot2_ns_crop.rds")
anc_LFC <- readRDS("results/parameter_sweeps/pLFC_S.pB10/pLFC_S.pB10_inv_change_plot2_ns_crop.rds")
mut_HFC <- readRDS("results/parameter_sweeps/pHFC_S.pB10-A/pHFC_S.pB10-A_inv_change_plot2_ns_crop.rds")
mut_LFC <- readRDS("results/parameter_sweeps/pLFC_S.pB10-A/pLFC_S.pB10-A_inv_change_plot2_ns_crop.rds")

# Crop figures ----
## Determine range ----
ph <- read_csv("input_data/strain_phenotypes.csv")

# Select strains to plot
ph_plot <- ph %>%
  filter(Strain %in% c("S.pB10", "S.pB10-A")) %>%
  mutate(log_Conj = log10(Conjugation_rate)) %>%
  mutate(Strain = case_when(Strain == "S.pB10" ~ anc_name,
                            Strain == "S.pB10-A" ~ mut_name))

anc_gamma <- ph_plot$log_Conj[ph_plot$Strain == anc_name]
mut_gamma <- ph_plot$log_Conj[ph_plot$Strain == mut_name]
anc_psi <- ph_plot$Growth_rate[ph_plot$Strain == anc_name]
mut_psi <- ph_plot$Growth_rate[ph_plot$Strain == mut_name]

mid_gamma <- mean(c(anc_gamma, mut_gamma))
mid_psi <- mean(c(anc_psi, mut_psi))

max_gamma <- mid_gamma + 3.5
min_gamma <- mid_gamma - 3.5

max_psi <- mid_psi + 0.35
min_psi <- mid_psi - 0.35



## Crop figures, remove legends and labels from inside plots
fig_theme <- theme_bw() +
  theme(legend.position = "top")

colors <- c(p_Anc, p_Mut)
shapes <- c(sh_Anc, sh_Mut)
names(colors) <- c(anc_name, mut_name)
names(shapes) <- c(anc_name, mut_name)

pA <- anc_LFC +
  geom_point(data = ph_plot,
             mapping = aes(log_Conj, Growth_rate, shape = Strain),
             size = psweep_point_size, color = point_color) +
  annotate("text",
           x = min_gamma + 0.4, y = max_psi - 0.04,
           label = "A", color = "white", size = label_size) +
  labs(title = paste0(anc_name, " as resident")) +
  scale_shape_manual(values = shapes) +
  scale_color_manual(values = colors) +
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust= 0.5)) +
  guides(shape = "none",
         color = "none")

pB <- mut_LFC +
  scale_y_continuous(limits = c(min_psi, max_psi),
                     expand = c(0.01, 0.01),
                     name = "LFC", position = "right") +
  geom_point(data = ph_plot,
             mapping = aes(log_Conj, Growth_rate, shape = Strain),
             size = psweep_point_size, color = point_color) +
  annotate("text",
           x = min_gamma + 0.4, y = max_psi - 0.04,
           label = "B", color = "white", size = label_size) +
  labs(title = paste0(mut_name, " as resident")) +
  scale_shape_manual(values = shapes) +
  scale_color_manual(values = colors) +
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(hjust= 0.5)) +
  guides(shape = "none",
         color = "none",
         fill = "none")

pC <- anc_HFC +
  geom_point(data = ph_plot,
             mapping = aes(log_Conj, Growth_rate, shape = Strain),
             size = psweep_point_size, color = point_color) +
  annotate("text",
           x = min_gamma + 0.4, y = max_psi - 0.04,
           label = "C", color = "black", size = label_size) +
  scale_shape_manual(values = shapes) +
  scale_color_manual(values = colors) +
  theme(legend.position = "none") +
  guides(shape = "none",
         color = "none",
         fill = "none")

pD <- mut_HFC +
  scale_y_continuous(limits = c(min_psi, max_psi),
                     expand = c(0.01, 0.01),
                     name = "HFC", position = "right") +
  geom_point(data = ph_plot,
             mapping = aes(log_Conj, Growth_rate, shape = Strain),
             size = psweep_point_size, color = point_color) +
  annotate("text",
           x = min_gamma + 0.4, y = max_psi - 0.04,
           label = "D", color = "black", size = label_size) +
  scale_shape_manual(values = shapes) +
  scale_color_manual(values = colors) +
  theme(axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.845, 0.73),
        legend.background = element_rect(fill = "white", color = "gray30",
                                         linewidth = 0.3),
        legend.title = element_text(size = 8)) +
  guides(shape = "none",
         color = "none",
         fill = guide_colourbar(display = "gradient",
                                ticks = FALSE,
                                theme = theme(legend.text = element_text(size = 5)),
                                barwidth = unit(0.15, "in"), barheight = unit(0.6, "in")))


final_plot <- (pA | pB) /
  (pC | pD) 

ggsave("figures/figS2_psweep_ns.pdf",
       final_plot, width = 5.5, height = 5.5, units = "in")



