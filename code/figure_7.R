source("code/functions.R")

a1 <- ggdraw() + draw_image("results/figures/d0_model_otu_Enterobacteriaceae (OTU 1).png")
a2 <- ggdraw() + draw_image("results/figures/d0_model_otu_Bacteroides (OTU 2).png")
a3 <- ggdraw() + draw_image("results/figures/d0_model_otu_Proteus (OTU 16).png")
a <- plot_grid(a1, a2, a3, labels = NULL, ncol = 1)
b1 <- ggdraw() + draw_image("results/figures/d0_model_otu_Lactobacillus (OTU 18).png")
b2 <- ggdraw() + draw_image("results/figures/d0_model_otu_Lachnospiraceae (OTU 9).png")
b3 <- ggdraw() + draw_image("results/figures/d0_model_otu_Clostridium (OTU 99).png")
b <- plot_grid(b1, b2, b3, labels = NULL, ncol = 1)
top_half <- plot_grid(a, b, labels = c("A", "B"), label_size = 12, ncol=2)
c <- ggdraw() + draw_image("results/figures/venn_overall_otus.png")

plot_grid(top_half, c, labels = c("", "C"), label_size = 12, ncol=1)+
  ggsave("results/figures/figure_7.pdf", width=6, height=8)+
  ggsave("submission/figure_7.pdf", width=6, height=8)

