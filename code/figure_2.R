source("code/functions.R")

a <- ggdraw() + draw_image("results/figures/richness_dn1.png")
c <- ggdraw() + draw_image("results/figures/richness_d0.png")
e <- ggdraw() + draw_image("results/figures/richness_d1.png")
b <- ggdraw() + draw_image("results/figures/shannon_dn1.png")
d <- ggdraw() + draw_image("results/figures/shannon_d0.png")
f <- ggdraw() + draw_image("results/figures/shannon_d1.png")


plot_grid(a, c, e, b, d, f, labels = c("A", "C", "E", "B", "D", "F"), label_size = 12, ncol=3)+
  ggsave("results/figures/figure_2.pdf", width=6, height=4)+
  ggsave("submission/figure_2.pdf", width=6, height=4)

