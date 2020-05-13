source("code/functions.R")

a <- ggdraw() + draw_image("results/figures/richness_dn1.png")
b <- ggdraw() + draw_image("results/figures/richness_d0.png")
c <- ggdraw() + draw_image("results/figures/richness_d1.png")
d <- ggdraw() + draw_image("results/figures/shannon_dn1.png")
e <- ggdraw() + draw_image("results/figures/shannon_d0.png")
f <- ggdraw() + draw_image("results/figures/shannon_d1.png")


plot_grid(a, b, c, d, e, f, labels = "AUTO", label_size = 12, ncol=3)+
  ggsave("results/figures/figure_3.pdf", width=6, height=4)+
  ggsave("submission/figure_3.pdf", width=6, height=4)

