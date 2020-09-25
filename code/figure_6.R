source("code/functions.R")

a <- ggdraw() + draw_image("results/figures/Bacteroides (OTU 2)_time.png")
b <- ggdraw() + draw_image("results/figures/Enterobacteriaceae (OTU 1)_time.png")
c <- ggdraw() + draw_image("results/figures/Enterococcus (OTU 23)_time.png")
d <- ggdraw() + draw_image("results/figures/Porphyromonadaceae (OTU 7)_time.png")
figure <- plot_grid(a, b, c, d, labels = "AUTO", label_size = 12, ncol=2)
legend <- ggdraw()+ draw_image("results/figures/overlap_OTUs_legend.png")


plot_grid(figure, legend, labels = NULL, label_size = 12, ncol=1, rel_heights = c(1, 0.2))+
  ggsave("results/figures/figure_6.pdf", width=6.875, height=4.3)+
  ggsave("submission/figure_6.tiff", width=6.875, height=4.3, dpi = 300, device = "tiff", compression = "lzw", units = "in")
