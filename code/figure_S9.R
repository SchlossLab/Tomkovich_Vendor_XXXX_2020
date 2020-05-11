source("code/functions.R")

a <- ggdraw() + draw_image("results/figures/Bacteroidaceae_time.png")
b <- ggdraw() + draw_image("results/figures/Deferribacteraceae.png")
c <- ggdraw() + draw_image("results/figures/Peptostreptococcaceae.png")
d <- ggdraw() + draw_image("results/figures/Porphyromonadaceae (OTU 7)_time.png")


plot_grid(a, b, c, d, labels = "AUTO", label_size = 12, ncol=1)+
  ggsave("results/figures/figure_S9.pdf", width=3, height=7)+
  ggsave("submission/figure_S9.pdf", width=3, height=7)
