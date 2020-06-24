source("code/functions.R")

a <- ggdraw() + draw_image("results/figures/Vendor_exp_scheme.png")
b <- ggdraw() + draw_image("results/figures/cfu_over_time.png")
c <- ggdraw() + draw_image("results/figures/weight_over_time.png")
d <- ggdraw() + draw_image("results/figures/C.diff_CFU_D7_stats.png")
e <- ggdraw() + draw_image("results/figures/weight_D2.png")
bottom_row <- plot_grid(d, e, labels = c("D", "E"), label_size = 12)

plot_grid(a, b, c, bottom_row, labels = c("A", "B", "C", ""), label_size = 12, ncol=1)+
  ggsave("results/figures/figure_1.pdf", width=3.875, height=8.5)+
  ggsave("submission/figure_1.pdf", width=3.875, height=8.5)


