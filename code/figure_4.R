source("code/functions.R")

a <- ggdraw() + draw_image("results/figures/richness_d1.png")
b <- ggdraw() + draw_image("results/figures/shannon_d1.png")
c <- ggdraw() + draw_image("results/figures/pcoa_day1.png")
e <- ggdraw() + draw_image("results/figures/dn1_vs_d7_thetayc.png")
left_panel <- plot_grid(a, b, c, e, labels = c("A", "B", "C", "E"), label_size = 12, ncol = 1)
d <- ggdraw() + draw_image("results/figures/D1top20_otus.png")

plot_grid(left_panel, d, labels = c("", "D"), label_size = 12, ncol=2, rel_heights = c(.9, 1.3), rel_widths = c(.9, 1.3))+
  ggsave("results/figures/figure_4.pdf", width=6, height=5.16)+
  ggsave("submission/figure_4.tiff", width=6, height=5.16, dpi = 600, device = "tiff", compression = "lzw", units = "in")
