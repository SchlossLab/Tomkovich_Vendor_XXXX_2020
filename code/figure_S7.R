source("code/functions.R")

a <- ggdraw() + draw_image("results/figures/Enterococcaceae_time.png")
b <- ggdraw() + draw_image("results/figures/Lachnospiraceae_time.png")


plot_grid(a, b, labels = "AUTO", label_size = 12, ncol=1)+
  ggsave("results/figures/figure_S7.pdf", width=3, height=3.5)+
  ggsave("submission/figure_S7.pdf", width=3, height=3.5)
