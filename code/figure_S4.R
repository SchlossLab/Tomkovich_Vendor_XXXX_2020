source("code/functions.R")

a <- ggdraw() + draw_image("results/figures/venn_dn1_otus.png")
b <- ggdraw() + draw_image("results/figures/venn_dn1_families.png")

plot_grid(a, b, labels = "AUTO", label_size = 12, ncol=1)+
  ggsave("results/figures/figure_S4.pdf", width=6, height=8)+
  ggsave("submission/figure_S4.pdf", width=6, height=8)
