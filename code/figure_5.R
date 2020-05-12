source("code/functions.R")

a <- ggdraw() + draw_image("results/figures/pcoa_day1.png")
b <- ggdraw() + draw_image("results/figures/Dn1toD1_otus_d1.png")
c <- ggdraw() + draw_image("results/figures/Dn1toD1_families_d1.png")


plot_grid(a, NULL, b, c, labels = c("A", "", "B", "C"), label_size = 12, ncol=2)+
  ggsave("results/figures/figure_5.pdf", width=6, height=6)+
  ggsave("submission/figure_5.pdf", width=6, height=6)
