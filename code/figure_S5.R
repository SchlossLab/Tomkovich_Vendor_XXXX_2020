source("code/functions.R")

a <- ggdraw() + draw_image("results/figures/venn_d0_otus.png")
b <- ggdraw() + draw_image("results/figures/venn_d0_families.png")

plot_grid(a, b, labels = "AUTO", label_size = 12, ncol=1)+
  ggsave("results/figures/figure_S5.pdf", width=6, height=8)+
  ggsave("submission/figure_S5.pdf", width=6, height=8)
