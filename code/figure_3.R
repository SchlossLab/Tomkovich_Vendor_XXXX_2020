source("code/functions.R")

a <- ggdraw() + draw_image("results/figures/pcoa_day0.png")
b <- ggdraw() + draw_image("results/figures/clind_impacted_families_plot.png")
c <- ggdraw() + draw_image("results/figures/richness_d0.png")
d <- ggdraw() + draw_image("results/figures/shannon_d0.png")
e <- ggdraw() + draw_image("results/figures/Dn1toD1_families_d0.png")
f <- ggdraw() + draw_image("results/figures/Dn1toD1_otus_d0.png")



plot_grid(a, b, c, d, e, f, labels = "AUTO", label_size = 12, ncol=2)+
  ggsave("results/figures/figure_3.pdf", width=6, height=8)+
  ggsave("submission/figure_3.pdf", width=6, height=8)

