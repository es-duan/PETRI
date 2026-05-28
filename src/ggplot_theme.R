# Common theme for ggplots

# Plot settings
fig_aes <- theme_bw() +
  theme(text = element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background =element_rect(fill="gray90"),
        strip.text.x = element_text(size = 12),
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10))

# Settings for parameter sweep plots
axes_aes <- theme_bw() +
  theme(text = element_text(size=20),
        axis.line = element_line(color = "black",linewidth = 1.5),
        axis.title = element_text(color = "black"),
        axis.text = element_text(size=1,color="black"),
        panel.border = element_blank())