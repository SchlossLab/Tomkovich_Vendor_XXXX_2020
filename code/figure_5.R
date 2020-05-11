source("code/functions.R")

a <- ggdraw() + draw_image("results/figures/class._logistic_regression_60.pdf")
b <- ggdraw() + draw_image("results/figures/class._logistic_regression_60_family.pdf")

plot_grid(a, b, labels = "AUTO", label_size = 12, ncol=1)+
  ggsave("results/figures/figure_5.pdf", width=6, height=8)+
  ggsave("submission/figure_5.pdf", width=6, height=8)
