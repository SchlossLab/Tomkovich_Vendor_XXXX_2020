source("code/functions.R")

a <- ggdraw() + draw_image("results/figures/clind_impacted_otus_plot_dn1.png")
b <- ggdraw() + draw_image("results/figures/clind_impacted_otus_plot_d0.png")
c <- ggdraw() + draw_image("results/figures/clind_impacted_families_plot_dn1.png")
d <- ggdraw() + draw_image("results/figures/clind_impacted_families_plot_d0.png")

plot_grid(a, b, c, d, labels = "AUTO", label_size = 12, ncol=2)+
  ggsave("results/figures/figure_5.pdf", width=8, height=11)+
  ggsave("submission/figure_5.pdf", width=8, height=11)

