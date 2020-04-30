source("code/functions.R")

a <- ggdraw() + draw_image("results/figures/exp1_cfu_time.png")
b <- ggdraw() + draw_image("results/figures/exp2_cfu_time.png")

plot_grid(a, b, labels = "AUTO", label_size = 12, ncol=1)+
  ggsave("results/figures/figure_S1.pdf", width=3.875, height=4.5)+
  ggsave("submission/figure_S1.pdf", width=3.875, height=4.5)


