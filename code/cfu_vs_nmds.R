source("code/functions.R")
library(vegan)

cfu <- metadata %>% select(id, vendor, day, count1, dilution1, count2, dilution2) %>%
  mutate(cfu1 = count1 / (20 ^ dilution1), cfu2 = count2 / (20 ^ dilution2)) %>%
  select(-starts_with("dilution"), -starts_with("count")) %>%
  gather(key = plate_num, value = cfu, -id, -vendor, -day) %>%
  filter(day == 1) %>%
  group_by(id) %>%
  summarize(mean_cfu = mean(cfu, na.rm = TRUE))

d1_data <- read_tsv("data/process/vendors.subsample.shared") %>%
  mutate(group = str_remove_all(Group, "No1")) %>%
  mutate(group = str_remove_all(Group, "No2")) %>%
  mutate(group = str_remove_all(Group, "no2")) %>%
  left_join(metadata, by = c("Group" = "id")) %>%
  filter(day == 1) %>%
  select(Group, starts_with("Otu")) %>%
  column_to_rownames(var = "Group")

d1_dist <- vegdist(d1_data, method = "bray")

d1_nmds <- metaMDS(d1_dist, k = 2) 

color_grad <- colorRampPalette(c("black", "red"))

d1_scores <- as.data.frame(scores(d1_nmds)) %>%
  mutate("group" = rownames(.)) %>%
  left_join(metadata, by = c("group" = "id")) %>%
  left_join(cfu, by = c("group" = "id")) %>%
  ggplot(aes(x = NMDS1, y = NMDS2, shape = vendor)) +
  geom_point(aes(color = mean_cfu)) +
  #scale_color_gradient(mean_cfu, low = "black", high = "red") +
  theme_light()
  
d1_scores
