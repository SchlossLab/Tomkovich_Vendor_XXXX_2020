source("code/functions.R")

#Drop any samples with NAs in the cfu column
cfu_data_final <- metadata %>% 
  filter(!is.na(cfu)) #378 observations that are not NAs

#Range of N mice per day
cfu_data_final %>% 
  group_by(day) %>% count() %>% arrange(n)

#Kruskal-wallis test of all timepoints where cfu data was collected (Days 0 through 9)
set.seed(19881117) #Same seed used for mothur analysis
cfu_kruskal_stats <- cfu_data_final %>% 
  select(day, vendor, cfu) %>% 
  group_by(day) %>% 
  nest() %>% 
  mutate(model=map(data, ~kruskal.test(x=.x$cfu, g=as.factor(.x$vendor)) %>%  tidy())) %>% 
  mutate(median = map(data, get_cfu_median_vendor)) %>% 
  unnest(c(model, median)) %>% 
  ungroup() #Ungroup before adjusting p values

#Adjust p-values for testing multiple days and write results to table
cfu_kruskal_stats_adjust <- cfu_kruskal_stats %>% 
  select(day, statistic, p.value, parameter, method, Schloss, Young, Jackson, `Charles River`, Taconic, Envigo) %>% 
  mutate(p.value.adj=p.adjust(p.value, method="BH")) %>% 
  arrange(p.value.adj) %>% 
  write_tsv("data/process/cfu_stats_all_days.tsv")

#List significant days after BH adjustment of p-values:
sig_cfu_days <- cfu_kruskal_stats_adjust %>% 
  filter(p.value.adj <= 0.05) %>% 
  pull(day)

#Perform pairwise Wilcoxan rank sum tests for days that were significant by Kruskal-Wallis rank sum test
cfu_stats_pairwise <- cfu_kruskal_stats %>% 
  filter(day %in% sig_cfu_days) %>% #only perform pairwise tests for days that were significant 
  group_by(day) %>% 
  mutate(model=map(data, ~pairwise.wilcox.test(x=.x$cfu, g=as.factor(.x$vendor), p.adjust.method="BH") %>% 
                     tidy() %>% 
                     mutate(compare=paste(group1, group2, sep="-")) %>% 
                     select(-group1, -group2) %>% 
                     pivot_wider(names_from=compare, values_from=p.value)
                   )
         ) %>% 
  unnest(model) %>% 
  select(-data, -parameter, -statistic) %>% 
  write_tsv("data/process/cfu_stats_sig_days.tsv")

#Format pairwise stats to use with ggpubr package
plot_format_stats <- cfu_stats_pairwise %>%
  #Remove all columns except pairwise comparisons and day
  select(-p.value, -method,-Schloss, -Young, -Jackson, -`Charles River`, -Taconic, -Envigo) %>% 
  group_split() %>% #Keeps a attr(,"ptype") to track prototype of the splits
  lapply(tidy_pairwise) %>% 
  bind_rows()

#Function to summarize data (calculate the median for each group) and plot the data
summarize_plot <- function(df){
  median_summary <- df %>% 
    group_by(vendor, day) %>% 
    summarize(median_cfu = median(cfu, na.rm = TRUE))
  ggplot(NULL) +
    geom_point(df, mapping = aes(x = day, y = cfu, color= vendor, fill = vendor, shape = experiment), size = 1.5, position = position_dodge(width = 0.6), show.legend = FALSE) +
    geom_line(median_summary, mapping = aes(x = day, y = median_cfu, color = vendor), size = 1.5) +
    scale_colour_manual(name=NULL,
                        values=color_scheme,
                        breaks=color_vendors,
                        labels=color_vendors)+
    scale_shape_manual(name="Experiment",
                       values=shape_scheme,
                       breaks=shape_experiment,
                       labels=shape_experiment) +
    labs(x = "Days Post-Infection", y = "CFU/g Feces") +
    scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9),
                       limits = c(-0.5, 9.5)) +
    geom_hline(yintercept = 100, linetype=2) +
    geom_text(x = 9.2, y = 102, color = "black", label = "LOD")+
    scale_y_log10(labels=fancy_scientific, breaks = c(10, 100, 10^3, 10^4, 10^5, 10^6, 10^7, 10^8, 10^9), limits = c(40, 350000000))+
    theme_classic()+
    theme(legend.position = "bottom",
          panel.grid.minor.x = element_line(size = 0.4, color = "grey"))  # Add gray lines to clearly separate symbols by days)
}

#CFU plot for the 1st experiment----
exp1_cfu <- summarize_plot(cfu_data_final %>% filter(experiment == 1))+
  ggtitle("1st experiment")+ #Title plot
  theme(plot.title = element_text(hjust = 0.5)) #Center plot title
save_plot(filename = paste0("results/figures/exp1_cfu_time.png"), plot = exp1_cfu, base_aspect_ratio = 2)

#CFU plot for the 2nd experiment----
exp2_cfu <- summarize_plot(cfu_data_final %>% filter(experiment == 2))+
  ggtitle("2nd experiment")+ #Title plot
  theme(plot.title = element_text(hjust = 0.5), #Center plot title
        legend.position = "none") 
save_plot(filename = paste0("results/figures/exp2_cfu_time.png"), plot = exp2_cfu, base_aspect_ratio = 2)

#CFU over time plot with astericks on days where cfu varied significantly across sources of mice using annotate()----
#List significant days after BH adjustment of p-values:
#Annotation labels determined by what days were significant
x_annotation <- sig_cfu_days 
y_position <- 10^8.5
label <- cfu_kruskal_stats_adjust %>% 
  mutate(p.value.adj=round(p.value.adj, digits = 4)) %>% 
  filter(p.value.adj <= 0.05) %>%
  mutate(p.signif = case_when(
    p.value.adj > 0.05 ~ "NS",
    p.value.adj <= 0.05 ~ "*"
  )) %>% 
  pull(p.signif)

#Plot data for all timepoints with annotated significance symbols for test of CFU differences across sources of mice at each timepoint
median_summary <- cfu_data_final %>% 
  group_by(vendor, day) %>% 
  summarize(median_cfu = median(cfu, na.rm = TRUE))
cfu_stats <- ggplot(NULL) +
  geom_point(cfu_data_final, mapping = aes(x = day, y = cfu, color= vendor, shape = experiment), size = 1.5, position = position_dodge(width = 0.6)) +
  geom_line(median_summary, mapping = aes(x = day, y = median_cfu, color = vendor), size = 1.5, show.legend = FALSE) +
  scale_colour_manual(name=NULL,
                      values=color_scheme,
                      breaks=color_vendors,
                      labels=color_vendors)+
  scale_shape_manual(name="Experiment",
                     values=shape_scheme,
                     breaks=shape_experiment,
                     labels=shape_experiment) +
  labs(x = "Days Post-Infection", y = "CFU/g Feces") +
  scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9),
                     limits = c(-0.5, 9.5)) +
  geom_hline(yintercept = 100, linetype=2) +
  geom_text(x = 9, y = 104, color = "black", label = "LOD")+
  scale_y_log10(labels=fancy_scientific, breaks = c(10, 100, 10^3, 10^4, 10^5, 10^6, 10^7, 10^8, 10^9), limits = c(NA, 10^9))+
  theme(text = element_text(size = 16))+  # Change font size for entire plot
  annotate("text", y = y_position, x = x_annotation, label = label, size =8)+
  theme_classic() +
  theme(legend.position = "bottom", #(0,0) bottom left (1,1) top right. #Previous: c(.92, .8)
        legend.key= element_rect(colour = "transparent", fill = "transparent"),
        panel.grid.minor.x = element_line(size = 0.4, color = "grey"))  # Add gray lines to clearly separate symbols by days)
save_plot(filename = paste0("results/figures/cfu_over_time.png"), plot = cfu_stats, base_aspect_ratio = 2)

#Plot of colonization status at day 7, separated by mouse colony sources
d7_title <- c(expression(paste(italic("C. difficile"), " status on Day 7"))) #Expression variable for the title so that bacteria name will be in italics
d7_status <- cfu_data_final %>% 
  filter(day == 7) %>% 
  mutate(clearance_status_d7 = case_when(clearance_status_d7 == "colonized" ~ "colonized", #Redo labeling so it makes sense on plot
                                         clearance_status_d7 == "not_detectable" ~ "cleared")) %>% 
  ggplot(aes(x=clearance_status_d7, group=vendor, fill=vendor))+
  scale_fill_manual(name=d7_title,
                      values=color_scheme,
                      breaks=color_vendors,
                      labels=color_vendors)+
  geom_bar()+
  facet_wrap(~vendor)+
  labs(title=d7_title,
       x=NULL,
       y="Number of mice")+ 
  theme_classic()+
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5)) # Center title
save_plot(filename = paste0("results/figures/c.diff_status_d7.png"), plot = d7_status, base_aspect_ratio = 2)

#Plot of colonization status at day 7, separated by day 7 colonization status
d7_title <- c(expression(paste(italic("C. difficile"), " status on Day 7"))) #Expression variable for the title so that bacteria name will be in italics
d7_status_binary <- cfu_data_final %>% 
  filter(day == 7) %>% 
  mutate(clearance_status_d7 = case_when(clearance_status_d7 == "colonized" ~ "colonized", #Redo labeling so it makes sense on plot
                                         clearance_status_d7 == "not_detectable" ~ "cleared")) %>% 
  ggplot(aes(x=clearance_status_d7, group=vendor))+
  geom_bar()+
  ylim(0, 25)+
  labs(title=d7_title,
       x=NULL,
       y="Number of mice")+ 
  theme_classic()+
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5)) # Center title
save_plot(filename = paste0("results/figures/c.diff_status_d7_binary.png"), plot = d7_status_binary, base_aspect_ratio = 2)

#Add clearance_status_dx columns for the rest of the days in addition to d7
metadata_cfu_status <- metadata %>% 
  mutate(clearance_status_d7 = case_when(cfu_d7 > 0 ~ "colonized",
                                         cfu_d7 == 0 ~ "cleared", #changed to match plot labeling used above
                                         cfu_d7 == NA ~ "NA"),
         clearance_status_d1 = case_when(cfu_d1 > 0 ~ "colonized",
                                         cfu_d1 == 0 ~ "cleared", 
                                         cfu_d1 == NA ~ "NA"),
         clearance_status_d2 = case_when(cfu_d2 > 0 ~ "colonized",
                                         cfu_d2 == 0 ~ "cleared", 
                                         cfu_d2 == NA ~ "NA"),
         clearance_status_d3 = case_when(cfu_d3 > 0 ~ "colonized",
                                         cfu_d3 == 0 ~ "cleared", 
                                         cfu_d3 == NA ~ "NA"),
         clearance_status_d4 = case_when(cfu_d4 > 0 ~ "colonized",
                                         cfu_d4 == 0 ~ "cleared", 
                                         cfu_d4 == NA ~ "NA"),
         clearance_status_d5 = case_when(cfu_d5 > 0 ~ "colonized",
                                         cfu_d5 == 0 ~ "cleared", 
                                         cfu_d5 == NA ~ "NA"),
         clearance_status_d6 = case_when(cfu_d6 > 0 ~ "colonized",
                                         cfu_d6 == 0 ~ "cleared", 
                                         cfu_d6 == NA ~ "NA"),
         clearance_status_d8 = case_when(cfu_d8 > 0 ~ "colonized",
                                         cfu_d8 == 0 ~ "cleared", 
                                         cfu_d8 == NA ~ "NA"),
         clearance_status_d9 = case_when(cfu_d9 > 0 ~ "colonized",
                                         cfu_d9 == 0 ~ "cleared", 
                                         cfu_d9 == NA ~ "NA"))


#Percent mice colonized each day, 2 experiments combined----
count_colonized <- function(timepoint, clearance_status_dx){
  metadata_cfu_status %>% 
    filter(day == timepoint, !is.na(cfu)) %>% #get rid of NA values
    count({{ clearance_status_dx }}) %>% 
    mutate(day = timepoint) %>% 
    rename(cfu_status = {{ clearance_status_dx }})
}  

#Figure out number of mice that either cleared C. difficile or were still colonized with C. difficile on each day
day_1 <- count_colonized(1, clearance_status_d1)
day_2 <- count_colonized(2, clearance_status_d2)
day_3 <- count_colonized(3, clearance_status_d3)
day_4 <- count_colonized(4, clearance_status_d4)
day_5 <- count_colonized(5, clearance_status_d5)
day_6 <- count_colonized(6, clearance_status_d6)
day_7 <- count_colonized(7, clearance_status_d7)
day_8 <- count_colonized(8, clearance_status_d8)
day_9 <- count_colonized(9, clearance_status_d9)

#Combine into a data frame and figure out percent colonized per day
percent_colonized_day <- rbind(day_1, day_2, day_3, day_4, day_5, day_6, day_7, day_8, day_9) %>% 
  group_by(day) %>% 
  mutate(total_n = sum(n)) %>% #Figure out total n per day (sum of cleared and colonized counts)
  filter(cfu_status == "colonized") %>% #Only need colonized rows to figure out percent colonized
  mutate(percent_colonized = (n/total_n) * 100)#Add new column for percent_colonized
#Day 6 and 7 are closest to 50% cleared and colonized. Day 7 has a larger N (40 vs 34)

#Plots of percent of mice colonized with C. difficile on each day for combined experiments and each experiment separately
# Function to plot % of mice colonized
plot_percent_colonized <- function(dataframe){
  dataframe %>% 
    ggplot(aes(x = day, y = percent_colonized)) +
    geom_col()+
    scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9),
                       limits = c(-0.5, 9.5)) +
    labs(title = NULL,
         x = "Days Post-infection", 
         y = "% of Mice Colonized")+
    theme_classic()+
    theme(text = element_text(size = 16)) + # Change font size for entire plot
    theme(plot.title = element_text(hjust = 0.5)) #Center plot title
  
}    
percent_colonized_plot <- plot_percent_colonized(percent_colonized_day)
save_plot(filename = paste0("results/figures/C.diff_percent_colonized.png"), percent_colonized_plot, base_height = 11, base_width = 8.5, base_aspect_ratio = 2)

# Plots of C. diff CFU data at timepoints where there are significant differences in CFU levels across the different sources of mice:

#Add p.value manually to timepoints of interest with ggpubr stat_pvalue_manual() function: Day 7 ----
#Figure out y.position via log transformation to match y-axis scale
y <- log10(max(cfu_data_final$cfu))

pairwise_wilcox_day7_plot <- plot_format_stats %>% 
  filter(day == 7) %>% 
  filter(p.adj <= 0.05) %>%  #Only show comparisons that were significant. p.value, which was adjusted < 0.05)
  mutate(p.adj="*") %>% #Just indicate whether statistically significant, exact p.adj values are in supplemental table
  mutate(y.position = c(y-1.0, y, y-0.5, y-1.4))

#Plot of Day 7 C. diff CFUs with p values added manually----
median_d7 <- cfu_data_final %>% 
  filter(day == 7) %>% 
  group_by(vendor) %>% 
  mutate(median_cfu = median(cfu, na.rm = TRUE)) %>% #create a column of median_cfu
  ungroup()
plot_CFU_D7 <- median_d7 %>% 
  ggplot(aes(x= vendor, y=cfu, color=vendor))+
  scale_colour_manual(name=NULL,
                      values=color_scheme,
                      breaks=color_vendors,
                      labels=color_vendors)+
  geom_hline(yintercept = 100, linetype=2) +
  geom_text(x = 11, y = 104, color = "black", label = "LOD")+
  geom_errorbar(aes(ymax = median_cfu, ymin = median_cfu), color = "gray50", size = 2)+ #Add lines to indicate the median for each group to the plot. Median calculated before y axis transformation
  geom_jitter(aes(shape=experiment, size=2)) +
  scale_shape_manual(name=NULL,
                     values=shape_scheme,
                     breaks=shape_experiment,
                     labels=shape_experiment) +
  labs(title="Day 7 Post-infection", 
       x=NULL,
       y="CFU/g Feces")+
  scale_y_log10(labels=fancy_scientific, breaks = c(10, 100, 10^3, 10^4, 10^5, 10^6, 10^7, 10^8, 10^9), limits = c(NA, 10^9))+
  scale_x_discrete(guide = guide_axis(n.dodge = 2))+  
  theme_classic()+
  theme(plot.title=element_text(hjust=0.5))+
  theme(legend.position = "none") + #Get rid of legend 
  stat_pvalue_manual(data = pairwise_wilcox_day7_plot, label = "p.adj", y.position = "y.position", size = 10, bracket.size = .8) +
  theme(text = element_text(size = 16)) # Change font size for entire plot
save_plot(filename = paste0("results/figures/C.diff_CFU_D7_stats.png"), plot_CFU_D7, base_height = 11, base_width = 8.5, base_aspect_ratio = 2)


