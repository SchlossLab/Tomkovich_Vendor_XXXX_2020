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
cfu_nas_final <- map(cfu_data_final, ~sum(is.na(.))) #217 for cfu1, 391 instances for cfu2. #182 for cfu

#Drop NAs from cfu column
cfu_data_final <- cfu_data_final %>% 
  filter(!is.na(cfu)) #381 observations that are not NAs

#Kruskal_wallis test for differences across groups at different timepoints with Benjamini-Hochburg correction. getOption("na.action") = "na.omit", so NAs are not included in statistical analysis----
kruskal_wallis_cfu <- cfu_data_final %>% 
  filter(day %in% c(1, 2, 3, 4, 5, 6, 7, 8, 9)) %>%  #only test days that we have CFU data for
  group_by(day) %>% 
  do(tidy(kruskal.test(cfu~factor(vendor), data=.))) %>% ungroup() %>% 
  mutate(p.value.adj=p.adjust(p.value, method="BH")) %>% 
  arrange(p.value.adj) 
#Timepoints where C. diff CFU is significantly different across the sources of mice
sig_C.diff_CFU_timepoints <- kruskal_wallis_cfu %>% 
    filter(p.value.adj <= 0.05) %>% 
    pull(day) 
#Days 5, 6, 7, 4, and 3 are when there are significant differences in C. difficile CFUs across the different sources of mice (listed in order of increasing adjusted P values)

#For significant timepoints, do pairwise.wilcox.test to determine which sources of mice are significantly different from each other regarding the amount of C. difficile CFUs.
pairwise.wilcox_groups <- function(timepoint){
  cfu_by_day <- cfu_data_final %>% 
    filter(day == timepoint)
  tidy(pairwise.wilcox.test(g = cfu_by_day$vendor, x = cfu_by_day$cfu, p.adjust.method = "BH"))
}

# Do pairwise.wilcox tests with BH correction for all significant timepoints----
for(d in sig_C.diff_CFU_timepoints){
  name <- paste("pairwise_wilcox_day", d, sep = "") #Way to name the data frames based on the date of interest
  assign(name, pairwise.wilcox_groups(d))
}
 pairwise_wilcox_day5
 pairwise_wilcox_day6
 pairwise_wilcox_day7
 pairwise_wilcox_day4
 pairwise_wilcox_day3

#Function to summarize data (calculate the mean for each group) and plot the data
summarize_plot <- function(df){
  mean_summary <- df %>% 
    group_by(vendor, day) %>% 
    summarize(mean_cfu = mean(cfu, na.rm = TRUE))
  ggplot(NULL) +
    geom_point(df, mapping = aes(x = day, y = cfu, color= vendor, fill = vendor), alpha = .2, size = .5, show.legend = FALSE, position = position_dodge(width = 0.6)) +
    geom_line(mean_summary, mapping = aes(x = day, y = mean_cfu, color = vendor), size = 1) +
    scale_colour_manual(name=NULL,
                        values=color_scheme,
                        breaks=color_vendors,
                        labels=color_vendors)+
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

# Boxplots of C. diff CFU data at timepoints where there are significant differences in CFU levels across the different sources of mice:
#Function to plot all significant genus relative abundances across vendors at a specific timepoint----
#Arguments:
# timepoint = timepoint to be analyzed
plot_C.diff_timepoint <- function(timepoint){
  plot_CFU_DX <- cfu_data_final %>% 
    filter(day == timepoint) %>% 
    ggplot(aes(x= vendor, y=cfu, color=vendor))+
    scale_colour_manual(name=NULL,
                        values=color_scheme,
                        breaks=color_vendors,
                        labels=color_vendors)+
    geom_hline(yintercept = 100, linetype=2) +
    geom_text(x = 11, y = 104, color = "black", label = "LOD")+
    geom_boxplot(outlier.shape = NA, size = 1.2)+
    geom_jitter(shape=19, size=2, alpha=0.6, position=position_jitterdodge(dodge.width=0.7, jitter.width=0.2)) +
    labs(title=NULL, 
         x=NULL,
         y="CFU/g Feces")+
    scale_y_log10(labels=fancy_scientific, breaks = c(10, 100, 10^3, 10^4, 10^5, 10^6, 10^7, 10^8, 10^9))+
    theme_classic()+
    theme(plot.title=element_text(hjust=0.5))+
    theme(legend.position = "none") + #Get rid of legend title & move legend position
    theme(text = element_text(size = 16))  #Remove legend
  save_plot(filename = paste0("results/figures/C.diff_CFU_D", timepoint,".png"), plot_CFU_DX, base_height = 11, base_width = 8.5, base_aspect_ratio = 2)
}
#Plot all the timepoints where C. diff CFUs were significantly different across sources of mice
for(d in sig_C.diff_CFU_timepoints){
  plot_C.diff_timepoint(d)
}

##Test of ggpubr----
library(ggpubr)
cfu_ggpubr <- cfu_data_final %>% 
  filter(day %in% c(1, 2, 3, 4, 5, 6, 7, 8, 9))  #only test days that we have CFU data for

KW_testcfu_ggpubr <- compare_means(cfu ~ vendor, data = cfu_ggpubr, method = kruskal.test, group.by = "day", p.adjust.method = "BH")
## Error
  
 
#Add p.value manually, also not working right
pairwise_wilcox_day5_plot <- pairwise_wilcox_day5 %>% 
  filter(p.value <= 0.05) %>%  #Only show comparisons that were significant. p.value, which was adjusted < 0.05)
  mutate(y.position = c(1320000, 1360000, 1400000, 1440000, 1480000, 1520000))
plot_CFU_D5 <- cfu_data_final %>% 
  filter(day == 5) %>% 
  ggplot(aes(x= vendor, y=cfu, color=vendor))+
  scale_colour_manual(name=NULL,
                      values=color_scheme,
                      breaks=color_vendors,
                      labels=color_vendors)+
  geom_hline(yintercept = 100, linetype=2) +
  geom_text(x = 11, y = 104, color = "black", label = "LOD")+
  geom_boxplot(outlier.shape = NA, size = 1.2)+
  geom_jitter(shape=19, size=2, alpha=0.6, position=position_jitterdodge(dodge.width=0.7, jitter.width=0.2)) +
  labs(title=NULL, 
       x=NULL,
       y="CFU/g Feces")+
  scale_y_log10(labels=fancy_scientific, breaks = c(10, 100, 10^3, 10^4, 10^5, 10^6, 10^7, 10^8, 10^9))+
  theme_classic()+
  theme(plot.title=element_text(hjust=0.5))+
  theme(legend.position = "none") + #Get rid of legend 
  theme(text = element_text(size = 16)) +  # Change font size for entire plot
  stat_pvalue_manual(data = pairwise_wilcox_day5_plot, label = "p.value", y.position = "y.position") 
save_plot(filename = paste0("results/figures/C.diff_CFU_D5_stats.png"), plot_CFU_D5, base_height = 11, base_width = 8.5, base_aspect_ratio = 2)
 ##Y-axis is off

#Add columns noting the C. diff CFU at significant timepoints for each mouse
cfu_data_final <-  cfu_data_final %>% ungroup()

cfu_d3 <- cfu_data_final %>% 
  filter(day == 3) %>% 
  mutate(cfu_d3 = cfu) %>% 
  select(mouse_id, cfu_d3)

cfu_d4 <- cfu_data_final %>% 
  filter(day == 4) %>% 
  mutate(cfu_d4 = cfu) %>% 
  select(mouse_id, cfu_d4)

cfu_d5 <- cfu_data_final %>% 
  filter(day == 5) %>% 
  mutate(cfu_d5 = cfu) %>% 
  select(mouse_id, cfu_d5)

cfu_d6 <- cfu_data_final %>% 
  filter(day == 6) %>% 
  mutate(cfu_d6 = cfu) %>% 
  select(mouse_id, cfu_d6)

cfu_d7 <- cfu_data_final %>% 
  filter(day == 7) %>% 
  mutate(cfu_d7 = cfu) %>% 
  select(mouse_id, cfu_d7)

#Merge all individual cfu_dx data frames onto cfu_data_final
cfu_data_final <- inner_join(cfu_data_final, cfu_d3, by = "mouse_id") 
cfu_data_final <- inner_join(cfu_data_final, cfu_d4, by = "mouse_id")
cfu_data_final <- inner_join(cfu_data_final, cfu_d5, by = "mouse_id")
cfu_data_final <- inner_join(cfu_data_final, cfu_d6, by = "mouse_id")
cfu_data_final <- inner_join(cfu_data_final, cfu_d7, by = "mouse_id")
 
