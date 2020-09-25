source("code/functions.R")

#Since A and B were faceted, need to add letter labels separately
letter_annotation <- data.frame(x_pos = c(.02, .55),
                                y_pos = c(.97, .97), 
                                lab = c("A", "B")) #To label the different panels of the plot

a <- ggdraw() + draw_image("results/figures/exp1_2_cfu_time.png")+
  geom_text(data = letter_annotation, aes( x = x_pos, y = y_pos, label = lab), size = 4.5, fontface = "bold")
c <- ggdraw() + draw_image("results/figures/C.diff_CFU_D7_stats.png")
d <- ggdraw() + draw_image("results/figures/weight_D2.png")
e <- ggdraw() + draw_image("results/figures/C.diff_percent_colonized.png")
bottom_panel <- plot_grid(c, d, e, labels = c("C", "D", "E"), label_size = 12, nrow = 1)



plot_grid(a, bottom_panel, labels = c("", ""), label_size = 12, ncol=1, rel_heights = c(1.4, 1))+
  ggsave("results/figures/figure_S2.pdf", width=6.875, height=5.16)+
  ggsave("submission/figure_S2.tiff", width=6.875, height=5.16, dpi = 300, device = "tiff", compression = "lzw", units = "in")


