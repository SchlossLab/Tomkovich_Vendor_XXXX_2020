source("code/functions.R")

a <- ggdraw() + draw_image("results/figures/c.diff_status_d7_binary.png")
b <- ggdraw() + draw_image("results/figures/d7_status_dn1_input.png")
c <- ggdraw() + draw_image("results/figures/d7_status_d0_input.png")
d <- ggdraw() + draw_image("results/figures/d7_status_d1_input.png")


plot_grid(a, b, c, d, labels = "AUTO", label_size = 12, ncol=2)+
  ggsave("results/figures/figure_S3.pdf", width=3.875, height=4.5)+
  ggsave("submission/figure_S3.pdf", width=3.875, height=4.5)
