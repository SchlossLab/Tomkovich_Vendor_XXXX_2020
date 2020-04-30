source("code/functions.R")

diversity_data <- read_tsv("data/process/vendors.diversity.summary") %>%
  mutate(group = str_remove_all(group, "No1")) %>%
  mutate(group = str_remove_all(group, "No2")) %>%
  mutate(group = str_remove_all(group, "no2")) %>%
  right_join(metadata, by = c("group" = "id"))

sobs_qq <- diversity_data %>%
  filter(!is.na(invsimpson), 
         method == "ave") %>%
  select(vendor, sobs) %>%
  ggplot(aes(sample = sobs, group = vendor, color = vendor)) +
  geom_qq() +
  stat_qq_line()

sobs_qq

invsimp_qq <- diversity_data %>%
  filter(!is.na(invsimpson), 
         method == "ave") %>%
  select(vendor, invsimpson) %>%
  mutate("log2_invsimp" = log2(invsimpson)) %>%
  ggplot(aes(sample = log2_invsimp, group = vendor, color = vendor)) +
  geom_qq() +
  stat_qq_line()

invsimp_qq

baseline_diversity <- diversity_data %>% 
  filter(!is.na(invsimpson), 
         method == "ave",
         day == -1) %>%
  select(vendor, sobs, invsimpson) %>%
  group_by(vendor)

sobs_boxplot <- ggplot(baseline_diversity,
                       aes(x = vendor, y = sobs,
                           color = vendor)) +
  geom_boxplot() +
  geom_jitter(shape = 19, size = 2) +
  scale_x_discrete(limits = c("Schloss", "Young", "Jackson",
                              "Charles River", "Taconic", "Envigo")) +
  labs(x = NULL, y = "Species Richness") +
  theme_classic()

sobs_boxplot
  
invsimp_boxplot <- baseline_diversity %>%
  mutate("log2_invsimp" = log2(invsimpson)) %>%
  ggplot(aes(x = vendor, y = log2_invsimp, color = vendor)) +
  geom_boxplot() +
  geom_jitter(shape = 19, size = 2) +
  scale_x_discrete(limits = c("Schloss", "Young", "Jackson",
                              "Charles River", "Taconic", "Envigo")) +
  labs(x = NULL, y = "log2 Inverse Simpson") +
  theme_classic()
  
invsimp_boxplot  

