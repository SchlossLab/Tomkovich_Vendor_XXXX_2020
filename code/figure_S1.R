source("code/functions.R")

a <- ggdraw() + draw_image("results/figures/pcoa_dn1_schloss.png")
b <- ggdraw() + draw_image("results/figures/pcoa_dn1_young.png")
c <- ggdraw() + draw_image("results/figures/pcoa_dn1_jackson.png")
d <- ggdraw() + draw_image("results/figures/pcoa_dn1_charles_river.png")
e <- ggdraw() + draw_image("results/figures/pcoa_dn1_taconic.png")
f <- ggdraw() + draw_image("results/figures/pcoa_dn1_envigo.png")
g <- ggdraw() + draw_image("results/figures/within_exp_dn1.png")
h <- ggdraw() + draw_image("results/figures/between_exp_dn1.png")

plot_grid(a, b, c, d, e, f, g, h, labels = "AUTO", label_size = 12, ncol=2)+
  ggsave("results/figures/figure_S1.pdf", width=5.4, height=9.0625)+
  ggsave("submission/figure_S1.tiff", width=5.4, height=9.0625, dpi = 300, device = "tiff", compression = "lzw", units = "in")



