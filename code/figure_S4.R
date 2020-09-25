source("code/functions.R")

a <- ggdraw() + draw_image("results/figures/venn_dn1_otus.png")
b <- ggdraw() + draw_image("results/figures/venn_d0_otus.png")
c <- ggdraw() + draw_image("results/figures/venn_d1_otus.png")


plot_grid(a, b, c, labels = "AUTO", label_size = 12, ncol=1)+
  ggsave("results/figures/figure_S4.pdf", width=6, height=8)+
  ggsave("submission/figure_S4.tiff", width=6, height=8, dpi = 300, device = "tiff", compression = "lzw", units = "in")
