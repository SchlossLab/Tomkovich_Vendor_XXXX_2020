source("code/functions.R")

a <- ggdraw() + draw_image("results/figures/Dn1toD1_otus_dn1.png")
b <- ggdraw() + draw_image("results/figures/Dn1toD1_otus_d0.png")
c <- ggdraw() + draw_image("results/figures/Dn1toD1_otus_d1.png")

plot_grid(a, b, c, labels = "AUTO", label_size = 12, ncol=3)+
  ggsave("results/figures/figure_5.pdf", width=8, height=3)+
  ggsave("submission/figure_5.pdf", width=8, height=3)

