source("code/functions.R")

cfu_plot <- metadata %>% select(vendor, day, count1, dilution1, count2, dilution2) %>%
  mutate(cfu1 = count1 / (20 ^ dilution1), cfu2 = count2 / (20 ^ dilution2)) %>%
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

#CFU quantification formula based on protocol used by the Young lab----
#Noticed 0 is counted as CFU even for dilutions where the -1 dilution was not checked. 
# If -1 dilution was not plated, we can not say for sure whether C. diff CFU for a mouse is really zero. Thus, zeros for any dilutions greater than the -1 dilution were converted to NAs
#Also need to be careful averaging counts together because if there's a 0 at a higher dilution (10^-6) but a count at a lower dilution (10^-5) the CFU count will be half of what it actually should be. (Above rule should take care of that).
# How to handle averaging...(hopefully rule to only keep 0s from -1 dilution will fix this issue)
cfu_data <- metadata %>% select(experiment, vendor, day, count1, dilution1, count2, dilution2, mouse_id, id) %>% 
  mutate(cfu1 = count1 * 20 * 1 / (10 ^ dilution1), cfu2 = count2 * 20 * (1 / (10 ^ dilution2))) %>% # Quantify CFU/g for each dilution that was plated
  select(-starts_with("count"))

#Need to figure out a way to make sure 0s are true 0s (i.e. only when the -1 was plated). 
#Number of 0s in cfu1 should equal number of -1s in dilution1. Since -1 dilution was never plated in dilution2 column, don't have to worry about counting any zeros from that column.
#Quantify how many instances we have 0s for cfu1
cfu1_0s <- cfu_data %>% filter(cfu1 == 0) #245/563 instances
cfu1_d <- cfu_data %>% filter(dilution1 == -1 & cfu1 == 0) #184 instances, meaning there should only be 184 0s.
cfu2_0s <- cfu_data %>% filter(cfu2 == 0) #133/563 instances. Most of these should be transformed to NAs
cfu_nas <- map(cfu_data, ~sum(is.na(.))) #156 for cfu1, 258 for cfu2

#Transform data so that 0s from -1 dilution remain 0s, but 0s for any dilutions beyond -1 become NA.----
cfu_data_final <- cfu_data %>%
  mutate(cfu1 = ifelse(cfu1 == 0 & dilution1 != -1, NA, cfu1)) %>% #Keeps 0s for -1 dilution, replaces 0s from any other dilution with NA
  mutate_at(vars(cfu2), ~replace(., . == 0, NA)) %>% #Changes all 0s for cfu2 to NA since the -1 dilution (the limit of detection) was not plated.
  group_by(id, experiment, mouse_id, day) %>% 
  mutate(cfu = mean(c(cfu1, cfu2), na.rm = TRUE)) %>% #Create a final cfu per ID (combination of mouse ID & date sample was collected) based on the average CFU/g for cfu1 and cfu2
  mutate(cfu = na_if(cfu, "NaN")) #Changes the NaNs in cfu column back to Nas
cfu1_0s <- cfu_data_final %>% filter(cfu1 == 0) #Now only 184 instances, which is what we predicted on line 21
cfu2_0s <- cfu_data_final %>% filter(cfu2 == 0) #0 instances of 0.
cfu_nas_final <- map(cfu_data_final, ~sum(is.na(.))) #217 for cfu1, 391 instances for cfu2

#Function to summarize data (calculate the mean for each group) and plot the data
summarize_plot <- function(df){
  mean_summary <- df %>% 
    group_by(vendor, day) %>% 
    summarize(mean_cfu = mean(cfu, na.rm = TRUE))
  ggplot(NULL) +
    geom_point(df, mapping = aes(x = day, y = cfu, color= vendor, fill = vendor), alpha = .2, size = .5, show.legend = FALSE, position = position_dodge(width = 0.6)) +
    geom_line(mean_summary, mapping = aes(x = day, y = mean_cfu, color = vendor), size = 1) +
    labs(x = "Days Post-Infection", y = "CFU/g Feces") +
    scale_x_continuous(breaks = c(-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
                       limits = c(-1.5, 10.5)) +
    geom_hline(yintercept = 100, linetype=2) +
    geom_text(x = 11, y = 104, color = "black", label = "LOD")+
    scale_y_log10(labels=fancy_scientific, breaks = c(10, 100, 10^3, 10^4, 10^5, 10^6, 10^7, 10^8, 10^9))+
    theme_classic()
}

#CFU plot that combines the 2 experiments----
combined_exp_cfu <- summarize_plot(cfu_data_final)

#CFU plot for the 1st experiment----
exp1_cfu <- summarize_plot(cfu_data_final %>% filter(experiment == 1))
#Note: weight and cfu data was only recorded for 1 mouse on D10 of experiment.

#CFU plot for the 2nd experiment----
exp2_cfu <- summarize_plot(cfu_data_final %>% filter(experiment == 2))

plot_grid(combined_exp_cfu, exp1_cfu, exp2_cfu, labels = c("Combined Experiments", "Experiment 1", "Experiment 2"), ncol = 1, label_x = .2, label_y = 1)+
  ggsave("exploratory/notebook/cfu_over_time.pdf", width = 8.5, height = 11)

