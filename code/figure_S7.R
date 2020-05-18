source("code/functions.R")

#Key families that vary with source and clindamycin treatment
a <- ggdraw() + draw_image("results/figures/Enterococcaceae_time.png")
b <- ggdraw() + draw_image("results/figures/Lachnospiraceae_time.png")
A <- plot_grid(a, b, labels = c("", ""), label_size = 12, nrow=1)

#Key families that vary with source
c <- ggdraw() + draw_image("results/figures/Bacteroidaceae_time.png")
d <- ggdraw() + draw_image("results/figures/Deferribacteraceae_time.png")
B <- plot_grid(c, d, labels = c("", ""), label_size = 12, nrow=1)

#Key families that vary with clindamycin treatment
e <- ggdraw() + draw_image("results/figures/Bifidobacteriaceae_time.png")
f <- ggdraw() + draw_image("results/figures/Coriobacteriaceae_time.png")
g <- ggdraw() + draw_image("results/figures/Ruminococcaceae_time.png")
h <- ggdraw() + draw_image("results/figures/Verrucomicrobiaceae_time.png")
C <- plot_grid(e, f, labels = c("", ""), label_size = 12, nrow=1)
D <- plot_grid(g, h, labels = c("", ""), label_size = 12, nrow=1)

#Combine all plots into one figure
plot_grid(A, B, C, D, labels = c("A", "B", "C", ""), label_size = 12, nrow=4, align = "t")+
  ggsave("results/figures/figure_S7.pdf", width=8, height=7)+
  ggsave("submission/figure_S7.pdf", width=8, height=7)
 
