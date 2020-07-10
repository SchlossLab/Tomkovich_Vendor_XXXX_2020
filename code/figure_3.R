source("code/functions.R")

a <- ggdraw() + draw_image("results/figures/richness_d0.png")
b <- ggdraw() + draw_image("results/figures/shannon_d0.png")
c <- ggdraw() + draw_image("results/figures/pcoa_day0.png")
top_row <- plot_grid(a, b, c, labels = "AUTO", label_size = 12, nrow=1)
d <- ggdraw() + draw_image("results/figures/clind_impacted_otus_plot.png")
e <- ggdraw() + draw_image("results/figures/D0top18_otus.png")
bottom_row <- plot_grid(d, e, labels = c("D", "E"), label_size = 12, nrow=1)

plot_grid(top_row, bottom_row, labels = c("", ""), label_size = 12, ncol=1)+
  ggsave("results/figures/figure_3.pdf", width=7.5, height=5.75)+
  ggsave("submission/figure_3.pdf", width=7.5, height=5.75)


