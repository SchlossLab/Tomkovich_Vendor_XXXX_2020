source("code/functions.R")

a <- ggdraw() + draw_image("results/figures/richness_d0.png")
b <- ggdraw() + draw_image("results/figures/shannon_d0.png")
c <- ggdraw() + draw_image("results/figures/pcoa_day0.png")
left_panel <- plot_grid(a, b, c, labels = "AUTO", label_size = 12, ncol=1)
d <- ggdraw() + draw_image("results/figures/D0top18_otus.png")
e <- ggdraw() + draw_image("results/figures/clind_impacted_otus_plot.png")
right_panel <- plot_grid(d, e, labels = c("D", "E"), label_size = 12, ncol=2, rel_widths = c(1, 1)) 
right_panel <- plot_grid(right_panel, NULL, labels = NULL, label_size = 12, nrow=2, rel_heights = c(1, 0.5)) #Add NULL plot to fix margin spacing in final figure

plot_grid(left_panel, right_panel, labels = c("", "D", "E"), label_size = 12, ncol=2, axis = "t", align = "v", rel_heights = c(.85, 3), rel_widths = c(.9, 2))+
  ggsave("results/figures/figure_3.pdf", width=6.87, height=5)+
  ggsave("submission/figure_3.tiff", width=6.87, height=5, dpi = 600, device = "tiff", compression = "lzw", units = "in")
  