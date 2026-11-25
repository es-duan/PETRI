# Export the 90% Dim simulation as Fig 7

# Load packages ----
library(tidyverse)

plot <- readRDS("figures/panels/Dim_psweep90.rds")

ggsave("figures/fig7_Dim_psweep.pdf",
       plot, width = 4, height = 3, units = "in")

