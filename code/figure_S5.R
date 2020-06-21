source("code/functions.R")

a <- ggdraw() + draw_image("results/figures/venn_dn1_families.png")
b <- ggdraw() + draw_image("results/figures/venn_d0_families.png")
c <- ggdraw() + draw_image("results/figures/venn_d1_families.png")
d <- ggdraw() + draw_image("results/figures/venn_overall_families.png")

plot_grid(a, b, c, d, labels = "AUTO", label_size = 12, ncol=1)+
  ggsave("results/figures/figure_S5.pdf", width=3, height=9.5)+
  ggsave("submission/figure_S5.pdf", width=3, height=9.5)
