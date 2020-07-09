source("code/functions.R")

a <- ggdraw() + draw_image("results/figures/Bacteroides (OTU 2)_time.png")
b <- ggdraw() + draw_image("results/figures/Enterobacteriaceae (OTU 1)_time.png")
c <- ggdraw() + draw_image("results/figures/Enterococcus (OTU 23)_time.png")
d <- ggdraw() + draw_image("results/figures/Porphyromonadaceae (OTU 7)_time.png")


plot_grid(a, b, c, d, labels = "AUTO", label_size = 12, ncol=2)+
  ggsave("results/figures/figure_8.pdf", width=7.5, height=4)+
  ggsave("submission/figure_8.pdf", width=7.5, height=4)
