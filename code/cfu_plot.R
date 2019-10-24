source("code/functions.R")

metadata <- read_csv("data/process/vendor_metadata.csv")

cfu_plot <- metadata %>% select(vendor, day, count1, dilution1, count2, dilution2) %>%
  mutate(cfu1 = count1 / (20 ^ dilution1), cfu2 = count2 / (20 ^ dilution2)) #%>%
  select(-starts_with("dilution"), -starts_with("count")) %>%
  gather(key = plate_num, value = cfu, -vendor, -day) %>%
  group_by(vendor, day) %>%
  summarize(mean_cfu = mean(cfu, na.rm = TRUE), sd_cfu = sd(cfu, na.rm = TRUE)) %>%
  mutate(se_cfu = calc_se(sd_cfu, n()),
         lower_ci = lower_ci(mean_cfu, se_cfu, n()),
         upper_ci = upper_ci(mean_cfu, se_cfu, n())) %>%
  ggplot(aes(x = day, y = mean_cfu, group = vendor, color = vendor)) +
  geom_line() +
  geom_point() +
  #geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci)) +
  scale_y_log10()

cfu_plot
