source("code/functions.R")

a1 <- ggdraw() + draw_image("results/figures/d0_model_otu_Enterobacteriaceae (OTU 1).png")
a2 <- ggdraw() + draw_image("results/figures/d0_model_otu_Bacteroides (OTU 2).png")
a3 <- ggdraw() + draw_image("results/figures/d0_model_otu_Proteus (OTU 16).png")
a <- plot_grid(a1, a2, a3, labels = NULL, ncol = 1)
top_half <- plot_grid(a, labels = c("A"), label_size = 12, ncol=1)
b <- ggdraw() + draw_image("results/figures/venn_overall_otus.png")

plot_grid(top_half, b, labels = c("", "B"), label_size = 12, ncol=1)+
  ggsave("results/figures/figure_5.pdf", width=6, height=9)+
  ggsave("submission/figure_5.pdf", width=6, height=9)

