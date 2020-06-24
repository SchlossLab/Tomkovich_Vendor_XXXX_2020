source("code/functions.R")

a1 <- ggdraw() + draw_image("results/figures/within_exp_dn1.png")
a2 <- ggdraw() + draw_image("results/figures/between_exp_dn1.png")
a3 <- ggdraw() + draw_image("results/figures/within_vendors_dn1.png")
a4 <- ggdraw() + draw_image("results/figures/between_vendors_dn1.png")
a <- plot_grid(a1, a2, a3, a4, labels = NULL, nrow = 1)
b1 <- ggdraw() + draw_image("results/figures/within_exp_d0.png")
b2 <- ggdraw() + draw_image("results/figures/between_exp_d0.png")
b3 <- ggdraw() + draw_image("results/figures/within_vendors_d0.png")
b4 <- ggdraw() + draw_image("results/figures/between_vendors_d0.png")
b <- plot_grid(b1, b2, b3, b4, labels = NULL, nrow = 1)
c1 <- ggdraw() + draw_image("results/figures/within_exp_d1.png")
c2 <- ggdraw() + draw_image("results/figures/between_exp_d1.png")
c3 <- ggdraw() + draw_image("results/figures/within_vendors_d1.png")
c4 <- ggdraw() + draw_image("results/figures/between_vendors_d1.png")
c <- plot_grid(c1, c2, c3, c4, labels = NULL, nrow = 1)

plot_grid(a, b, c, labels = c("A", "B", "C"), label_size = 12, ncol=1)+
  ggsave("results/figures/figure_4.pdf", width=6, height=5)+
  ggsave("submission/figure_4.pdf", width=6, height=5)
