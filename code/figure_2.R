source("code/functions.R")

a <- ggdraw() + draw_image("results/figures/pcoa_day-1.png")
b <- ggdraw() + draw_image("results/figures/richness_dn1.png")
c <- ggdraw() + draw_image("results/figures/shannon_dn1.png")
d <- ggdraw() + draw_image("results/figures/Dn1toD1_families_dn1.png")
e <- ggdraw() + draw_image("results/figures/Dn1toD1_otus_dn1.png")



plot_grid(a, NULL, b, c, d, e,  labels = c("A", "", "B", "C", "D", "E"), label_size = 12, ncol=2)+
  ggsave("results/figures/figure_2.pdf", width=6, height=8)+
  ggsave("submission/figure_2.pdf", width=6, height=8)


