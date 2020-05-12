source("code/functions.R")

a <- ggdraw() + draw_image("results/figures/clind_impacted_otus_plot.png")
b <- ggdraw() + draw_image("results/figures/clind_impacted_families_plot.png")

plot_grid(a, b, labels = "AUTO", label_size = 12, ncol=2)+
  ggsave("results/figures/figure_4.pdf", width=7, height=6.2)+
  ggsave("submission/figure_4.pdf", width=7, height=6.2)

