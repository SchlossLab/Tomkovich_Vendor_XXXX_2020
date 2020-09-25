source("code/functions.R")

a <- ggdraw() + draw_image("results/figures/Vendor_exp_scheme.png")
b <- ggdraw() + draw_image("results/figures/cfu_over_time.png")
c <- ggdraw() + draw_image("results/figures/weight_over_time.png")


plot_grid(a, b, c, labels = c("A", "B", "C", ""), label_size = 12, ncol=1)+
  ggsave("results/figures/figure_2.pdf", width=5, height=7.5)+
  ggsave("submission/figure_2.tiff", width=5, height=7.5, dpi = 300, device = "tiff", compression = "lzw", units = "in")



