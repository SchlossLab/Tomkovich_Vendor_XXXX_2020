source("code/functions.R")

weight_plot <- metadata %>% select(vendor, day, weight) %>%
  group_by(vendor, day) %>%
  summarize(mean_weight = mean(weight, na.rm = TRUE), 
            sd_weight = sd(weight, na.rm = TRUE)) %>%
  mutate(se_weight = calc_se(sd_weight, n()),
         lower_ci = lower_ci(mean_weight, se_weight, n()),
         upper_ci = upper_ci(mean_weight, se_weight, n())) %>%
  ggplot(aes(x = day, y = mean_weight, group = vendor, color = vendor)) +
  geom_line() +
  geom_errorbar(aes(ymax = upper_ci, ymin = lower_ci, width = 0.5)) +
  geom_point() +
  labs(x = "Day", y = "Mean Weight (g)") +
  scale_x_continuous(breaks = c(-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
                     limits = c(-1.5, 10.5)) +
  coord_fixed() +
  theme_light()

weight_plot
