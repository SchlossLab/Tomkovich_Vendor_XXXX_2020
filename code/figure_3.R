source("code/functions.R")

a <- ggdraw() + draw_image("results/figures/pcoa_day-1.png")
b <- ggdraw() + draw_image("results/figures/pcoa_day0.png")
c <- ggdraw() + draw_image("results/figures/pcoa_day1.png")

plot_grid(a, b, c, labels = "AUTO", label_size = 12, ncol=1)+
  ggsave("results/figures/figure_3.pdf", width=3.4, height=8.5)+
  ggsave("submission/figure_3.pdf", width=3.4, height=8.5)


