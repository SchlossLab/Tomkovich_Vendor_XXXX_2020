source("code/functions.R")

a <- ggdraw() + draw_image("results/figures/exp1_cfu_time.png")
b <- ggdraw() + draw_image("results/figures/exp2_cfu_time.png")
c <- ggdraw() + draw_image("results/figures/C.diff_percent_colonized.png")

plot_grid(a, b, c, labels = "AUTO", label_size = 12, ncol=2)+
  ggsave("results/figures/figure_S1.pdf", width=8, height=5.5)+
  ggsave("submission/figure_S1.pdf", width=8, height=5.5)


