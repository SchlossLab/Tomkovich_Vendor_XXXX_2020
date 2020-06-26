source("code/functions.R")

a <- ggdraw() + draw_image("results/figures/pcoa_dn1_schloss.png")
b <- ggdraw() + draw_image("results/figures/pcoa_dn1_young.png")
c <- ggdraw() + draw_image("results/figures/pcoa_dn1_jackson.png")
d <- ggdraw() + draw_image("results/figures/pcoa_dn1_charles_river.png")
e <- ggdraw() + draw_image("results/figures/pcoa_dn1_taconic.png")
f <- ggdraw() + draw_image("results/figures/pcoa_dn1_envigo.png")

plot_grid(a, b, c, d, e, f, labels = "AUTO", label_size = 12, ncol=2)+
  ggsave("results/figures/figure_S2.pdf", width=6, height=8)+
  ggsave("submission/figure_S2.pdf", width=6, height=8)


