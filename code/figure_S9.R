source("code/functions.R")

a <- ggdraw() + draw_image("results/figures/Bifidobacteriaceae_time.png")
b <- ggdraw() + draw_image("results/figures/Coriobacteriaceae_time.png")
c <- ggdraw() + draw_image("results/figures/Ruminococcaceae_time.png")
d <- ggdraw() + draw_image("results/figures/Verrucomicrobiaceae_time.png")


plot_grid(a, b, c, d, labels = "AUTO", label_size = 12, ncol=1)+
  ggsave("results/figures/figure_S9.pdf", width=3, height=7)+
  ggsave("submission/figure_S9.pdf", width=3, height=7)
