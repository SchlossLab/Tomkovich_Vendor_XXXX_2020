source("code/functions.R")

a <- ggdraw() + draw_image("results/figures/richness_dn1.png")
b <- ggdraw() + draw_image("results/figures/shannon_dn1.png")
c <- ggdraw() + draw_image("results/figures/pcoa_day-1.png")
left_panel <- plot_grid(a, b, c, labels = c("A", "B", "C"), label_size = 12, ncol = 1)
d <- ggdraw() + draw_image("results/figures/Dn1top20_otus.png")

plot_grid(left_panel, d, labels = c("", "D"), label_size = 12, ncol=2, rel_heights = c(.9, 1.2), rel_widths = c(.9, 1.2))+
  ggsave("results/figures/figure_1.pdf", width=6.8, height=4.73)+
  ggsave("submission/figure_1.tiff", width=6.8, height=4.73, dpi = 600, device = "tiff", compression = "lzw", units = "in")


