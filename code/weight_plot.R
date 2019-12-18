source("code/functions.R")

#Combined experiments----
weight_data <- metadata %>% select(experiment, vendor, day, weight)

mean_weight_summary <- weight_data %>% 
  group_by(vendor, day) %>%
  summarize(mean_weight = mean(weight, na.rm = TRUE))

weight_plot <- ggplot(NULL) +
  geom_point(weight_data, mapping = aes(x = day, y = weight, color= vendor, fill = vendor), alpha = .2, size = .5, show.legend = FALSE, position = position_dodge(width = 0.6)) +
  geom_line(mean_weight_summary, mapping = aes(x = day, y = mean_weight, color = vendor), size = 1) +
  labs(x = "Day", y = "Weight (g)") +
  scale_x_continuous(breaks = c(-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
                     limits = c(-1.5, 10.5)) +
  coord_fixed() +
  theme_classic()

#Experiment 1 ----
weight_data <- metadata %>% select(experiment, vendor, day, weight) %>% 
  filter(experiment == 1)

mean_weight_summary <- weight_data %>% 
  group_by(vendor, day) %>%
  summarize(mean_weight = mean(weight, na.rm = TRUE))

weight_plot_exp_1 <- ggplot(NULL) +
  geom_point(weight_data, mapping = aes(x = day, y = weight, color= vendor, fill = vendor), alpha = .2, size = .5, show.legend = FALSE, position = position_dodge(width = 0.6)) +
  geom_line(mean_weight_summary, mapping = aes(x = day, y = mean_weight, color = vendor), size = 1) +
  labs(x = "Day", y = "Weight (g)") +
  scale_x_continuous(breaks = c(-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
                     limits = c(-1.5, 10.5)) +
  coord_fixed() +
  theme_classic()

#Experiment 2 ----
weight_data <- metadata %>% select(experiment, vendor, day, weight) %>% 
  filter(experiment == 2)

mean_weight_summary <- weight_data %>% 
  group_by(vendor, day) %>%
  summarize(mean_weight = mean(weight, na.rm = TRUE))

weight_plot_exp_2 <- ggplot(NULL) +
  geom_point(weight_data, mapping = aes(x = day, y = weight, color= vendor, fill = vendor), alpha = .2, size = .5, show.legend = FALSE, position = position_dodge(width = 0.6)) +
  geom_line(mean_weight_summary, mapping = aes(x = day, y = mean_weight, color = vendor), size = 1) +
  labs(x = "Day", y = "Weight (g)") +
  scale_x_continuous(breaks = c(-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
                     limits = c(-1.5, 10.5)) +
  coord_fixed() +
  theme_classic()

#Graph weight data as percent baseline weight (see Theriot et al. 2011).----
#Percent baseline weight for each mouse will be calculated for each mouse
#based on the weight recorded on D-1 of the experiment

#Create data frame that has just the baseline_weight data for D-1 of each mouse
baseline_weight <- metadata %>% select(mouse_id, weight, day) %>% 
  filter(day == -1) %>% 
  mutate(baseline_weight = weight) %>% 
  select(mouse_id, baseline_weight)

#Join baseline data frame to main metadata
baseline_weight_data <- inner_join(metadata, baseline_weight, by = "mouse_id") %>% 
  select(experiment, mouse_id, vendor, day, weight, baseline_weight)

#Calculate percent baseline weight data for each mouse based on the D-1 weight.
percent_baseline_weight_data <- baseline_weight_data %>%  
  group_by(mouse_id, day) %>% 
  mutate(percent_baseline_weight = 100 + ((weight-baseline_weight)/(baseline_weight))*100)

#Function to summarize data (calculate the mean for each group) and plot the data
summarize_plot <- function(df){
  mean_summary <- df %>% 
    group_by(vendor, day) %>% 
    summarize(mean_percent_weight = mean(percent_baseline_weight, na.rm = TRUE))
  ggplot(NULL) +
    geom_point(df, mapping = aes(x = day, y = percent_baseline_weight, color= vendor, fill = vendor), alpha = .2, size = .5, show.legend = FALSE, position = position_dodge(width = 0.6)) +
    geom_line(mean_summary, mapping = aes(x = day, y = mean_percent_weight, color = vendor), size = 1) +
    labs(x = "Days Post-Infection", y = "% Baseline Weight") +
    scale_x_continuous(breaks = c(-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
                       limits = c(-1.5, 10.5)) +
    theme_classic()
}

#Combined percent baseline weight plot for the 2 experiments----
combined_exp_weight <- summarize_plot(percent_baseline_weight_data)

#Percent baseline weight plot for the 1st experiment----
exp1_weight <- summarize_plot(percent_baseline_weight_data %>% filter(experiment == 1))
#Note: weight and cfu data was only recorded for 1 mouse on D10 of experiment.

#Percent baseline weight plot for the 2nd experiment----
exp2_weight <- summarize_plot(percent_baseline_weight_data %>% filter(experiment == 2))

plot_grid(combined_exp_weight, exp1_weight, exp2_weight, labels = c("Combined Experiments", "Experiment 1", "Experiment 2"), ncol = 1, label_x = 0, label_y = 1)+
  ggsave("exploratory/notebook/weight_changes.pdf", width = 8.5, height = 11)


