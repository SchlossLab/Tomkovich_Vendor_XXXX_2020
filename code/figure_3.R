source("code/functions.R")

a <- ggdraw() + draw_image("results/figures/richness_d0.png")
b <- ggdraw() + draw_image("results/figures/shannon_d0.png")
c <- ggdraw() + draw_image("results/figures/pcoa_day0.png")
left_panel <- plot_grid(a, b, c, labels = "AUTO", label_size = 12, ncol=1)
d <- ggdraw() + draw_image("results/figures/D0top18_otus.png")
e <- ggdraw() + draw_image("results/figures/clind_impacted_otus_plot.png")

plot_grid(left_panel, d, e, labels = c("", "D", "E"), label_size = 12, ncol=3)+
  ggsave("results/figures/figure_3.pdf", width=8, height=4.25)+
  ggsave("submission/figure_3.pdf", width=8, height=4.25)


