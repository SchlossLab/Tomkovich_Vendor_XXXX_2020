source("code/functions.R")

a <- ggdraw() + draw_image("results/figures/c.diff_status_d7_binary.png")
b <- ggdraw() + draw_image("results/figures/c.diff_status_d7.png")
c <- ggdraw() + draw_image("results/figures/class._logistic_regression_60.png")

plot_grid(a, b, c, labels = "AUTO", label_size = 12, ncol=2)+
  ggsave("results/figures/figure_S3.pdf", width=6.875, height=5.04)+
  ggsave("submission/figure_S3.tiff", width=6.875, height=5.04, dpi = 300, device = "tiff", compression = "lzw", units = "in")
 
