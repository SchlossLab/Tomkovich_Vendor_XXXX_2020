source("code/functions.R")

a <- ggdraw() + draw_image("results/figures/exp1_cfu_time.png")
b <- ggdraw() + draw_image("results/figures/exp2_cfu_time.png")
c <- ggdraw() + draw_image("results/figures/C.diff_CFU_D7_stats.png")
d <- ggdraw() + draw_image("results/figures/weight_D2.png")
e <- ggdraw() + draw_image("results/figures/C.diff_percent_colonized.png")
right_panel <- plot_grid(c, d, e, labels = c("C", "D", "E"), label_size = 12, nrow = 1)


plot_grid(a, right_panel, b,  labels = c("A", "", "B"), label_size = 12, ncol=2)+
  ggsave("results/figures/figure_S1.pdf", width=8, height=4.5)+
  ggsave("submission/figure_S1.pdf", width=8, height=4.5)


