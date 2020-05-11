source("code/functions.R")

a <- ggdraw() + draw_image("results/figures/Bacteroidaceae_time.png")
b <- ggdraw() + draw_image("results/figures/Deferribacteraceae_time.png")
c <- ggdraw() + draw_image("results/figures/Peptostreptococcaceae_time.png")


plot_grid(a, b, c, labels = "AUTO", label_size = 12, ncol=1)+
  ggsave("results/figures/figure_S8.pdf", width=3, height=5)+
  ggsave("submission/figure_S8.pdf", width=3, height=5)
 