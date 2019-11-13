source("code/functions.R")
library(vegan)
library(cowplot)

dn1_data <- read_tsv("data/process/vendors.subsample.shared") %>%
  mutate(group = str_remove_all(Group, "No1")) %>%
  mutate(group = str_remove_all(Group, "No2")) %>%
  mutate(group = str_remove_all(Group, "no2")) %>%
  left_join(metadata, by = c("Group" = "id"))

exp1_data <- dn1_data %>%
  filter(day == -1, experiment == 1) %>%
  select(Group, starts_with("Otu")) %>%
  column_to_rownames(var = "Group")

exp2_data <- dn1_data %>%
  filter(day == -1, experiment == 2) %>%
  select(Group, starts_with("Otu")) %>%
  column_to_rownames(var = "Group")

exp1_dist <- vegdist(exp1_data, method = "bray")
exp2_dist <- vegdist(exp2_data, method = "bray")

exp1_nmds <- metaMDS(exp1_dist, k = 2) 
exp2_nmds <- metaMDS(exp2_dist, k = 2)

exp1_scores <- as.data.frame(scores(exp1_nmds)) %>%
  mutate("group" = rownames(.)) %>%
  left_join(metadata, by = c("group" = "id")) %>%
  ggplot(aes(x = NMDS1, y = NMDS2, color = vendor)) +
  geom_point() +
  theme_light()
  
exp2_scores <- as.data.frame(scores(exp2_nmds)) %>%
  mutate("group" = rownames(.)) %>%
  left_join(metadata, by = c("group" = "id")) %>%
  ggplot(aes(x = NMDS1, y = NMDS2, color = vendor)) +
  geom_point() +
  theme_light()

plot_grid(exp1_scores, exp2_scores, labels = c("A", "B"))
