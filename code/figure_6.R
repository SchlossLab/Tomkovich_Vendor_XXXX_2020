source("code/functions.R")

a <- ggdraw() + draw_image("results/figures/clind_impacted_otus_plot_dn1.png")
b <- ggdraw() + draw_image("results/figures/clind_impacted_otus_plot_d0.png")

plot_grid(a, b, labels = "AUTO", label_size = 12, ncol=2)+
  ggsave("results/figures/figure_6.pdf", width=8, height=5.5)+
  ggsave("submission/figure_6.pdf", width=8, height=5.5)
