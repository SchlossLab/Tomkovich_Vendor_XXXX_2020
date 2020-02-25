source("code/functions.R")

#Drop any samples with NAs in the cfu column
cfu_data_final <- metadata %>% 
  filter(!is.na(cfu)) #378 observations that are not NAs

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
 pairwise_wilcox_day5 %>% filter(p.value <= 0.05) # 6 sig. pairwise: Jax vs Schloss, Tac. vs Schloss, Jax vs Young, Tac. vs Young, CR vs Jax, CR vs Tac
 pairwise_wilcox_day6 %>% filter(p.value <= 0.05) # 6 sig. pairwise:
 pairwise_wilcox_day7 %>% filter(p.value <= 0.05) # 4 sig. pairwise:
 pairwise_wilcox_day4 %>% filter(p.value <= 0.05) # 0 sig. pairwise:
 pairwise_wilcox_day3 %>% filter(p.value <= 0.05) # 3 sig. pairwise:

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
    scale_x_continuous(breaks = c(-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9),
                       limits = c(-1.5, 9.5)) +
    geom_hline(yintercept = 100, linetype=2) +
    geom_text(x = 11, y = 104, color = "black", label = "LOD")+
    scale_y_log10(labels=fancy_scientific, breaks = c(10, 100, 10^3, 10^4, 10^5, 10^6, 10^7, 10^8, 10^9))+
    theme_classic()
}

#CFU plot that combines the 2 experiments----
combined_exp_cfu <- summarize_plot(cfu_data_final)

#CFU plot for the 1st experiment----
exp1_cfu <- summarize_plot(cfu_data_final %>% filter(experiment == 1))

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
    geom_jitter(shape=19, size=2, alpha=0.6, position=position_jitterdodge(dodge.width=0.5, jitter.width=0.2)) +
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
 
