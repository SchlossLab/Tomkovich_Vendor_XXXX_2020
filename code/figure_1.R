source("code/functions.R")

a <- ggdraw() + draw_image("results/figures/Vendor_exp_scheme.png")
b <- ggdraw() + draw_image("results/figures/cfu_over_time.png")
c <- ggdraw() + draw_image("results/figures/weight_over_time.png")
d <- ggdraw() + draw_image("results/figures/c.diff_status_d7.png")

plot_grid(a, b, c, d, labels = "AUTO", label_size = 12, ncol=1)+
  ggsave("results/figures/figure_1.pdf", width=3.875, height=8.5)+
  ggsave("submission/figure_1.pdf", width=3.875, height=8.5)


