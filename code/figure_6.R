source("code/functions.R")

a <- ggdraw() + draw_image("results/figures/clind_impacted_otus_plot.png")

plot_grid(a, labels = "", label_size = 12, ncol=1)+
  ggsave("results/figures/figure_6.pdf", width=8, height=5.5)+
  ggsave("submission/figure_6.pdf", width=8, height=5.5)
