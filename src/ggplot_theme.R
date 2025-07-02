# Common theme for ggplots

# Plot settings
fig_aes <- theme_bw() +
  theme(text = element_text(size=24),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background =element_rect(fill="gray90"),
        strip.text.x = element_text(size = 24),
        axis.title.x = element_text(size=30),
        axis.title.y = element_text(size=30),
        axis.text.x = element_text(size=24),
        axis.text.y = element_text(size=24))
