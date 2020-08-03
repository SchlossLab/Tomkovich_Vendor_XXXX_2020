source("code/functions.R")

a <- ggdraw() + draw_image("results/figures/d0_top10_otus.png")
b <- ggdraw() + draw_image("results/figures/venn_overall_otus.png")

plot_grid(a, b, labels = "AUTO", label_size = 12, ncol=1)+
  ggsave("results/figures/figure_5.pdf", width=7, height=9)+
  ggsave("submission/figure_5.pdf", width=7, height=9)

