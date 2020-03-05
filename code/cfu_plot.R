source("code/functions.R")

#Drop any samples with NAs in the cfu column
cfu_data_final <- metadata %>% 
  filter(!is.na(cfu)) #378 observations that are not NAs

#Kruskal-wallis test of all timepoints where cfu data was collected (Days 0 through 9)
set.seed(19881117) #Same seed used for mothur analysis
cfu_kruskal_stats <- cfu_data_final %>% 
  select(day, vendor, cfu) %>% 
  group_by(day) %>% 
  nest() %>% 
  mutate(model=map(data, ~kruskal.test(x=.x$cfu, g=as.factor(.x$vendor)) %>%  tidy())) %>% 
  mutate(mean = map(data, get_cfu_mean_vendor)) %>% 
  unnest(c(model, mean)) %>% 
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

#Function to tidy pairwise comparisons to use for adding stats to plots
tidy_pairwise <- function(spread_pairwise){
  spread_pairwise %>% 
    pivot_longer(-day, names_to = "compare", values_to = "p.adj") %>% 
    separate(col = compare, c("group1", "group2"), sep = "-", remove = TRUE)
}

#Format pairwise stats to use with ggpubr package
plot_format_stats <- cfu_stats_pairwise %>%
  #Remove all columns except pairwise comparisons and day
  select(-p.value, -method,-Schloss, -Young, -Jackson, -`Charles River`, -Taconic, -Envigo) %>% 
  group_split() %>% #Keeps a attr(,"ptype") to track prototype of the splits
  lapply(tidy_pairwise) %>% 
  bind_rows()

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
    scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9),
                       limits = c(-0.5, 9.5)) +
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
    scale_y_log10(labels=fancy_scientific, breaks = c(10, 100, 10^3, 10^4, 10^5, 10^6, 10^7, 10^8, 10^9), limits=c(NA, 10^9))+
    theme_classic()+
    theme(plot.title=element_text(hjust=0.5))+
    theme(legend.position = "none") + #Get rid of legend title & move legend position
    theme(text = element_text(size = 16))  #Remove legend
  save_plot(filename = paste0("results/figures/C.diff_CFU_D", timepoint,".png"), plot_CFU_DX, base_height = 11, base_width = 8.5, base_aspect_ratio = 2)
}
#Plot all the timepoints where C. diff CFUs were significantly different across sources of mice
for(d in sig_cfu_days){
  plot_C.diff_timepoint(d)
}

#Add p.value manually to timepoints of interest with ggpubr stat_pvalue_manual() function: Day 5 & 7 ----
#Figure out y.position via log transformation to match y-axis scale
y <- log10(max(cfu_data_final$cfu))
#Data frames of p.values to add manually
pairwise_wilcox_day5_plot <- plot_format_stats %>% 
  filter(day == 5) %>% 
  filter(p.adj <= 0.05) %>%  #Only show comparisons that were significant. p.value, which was adjusted < 0.05)
  mutate(p.adj=round(p.adj, digits = 4)) %>% 
  mutate(y.position = c(y-2, y, y-1.4, y-0.5, y-1.0, y-1.7))

pairwise_wilcox_day7_plot <- plot_format_stats %>% 
  filter(day == 7) %>% 
  filter(p.adj <= 0.05) %>%  #Only show comparisons that were significant. p.value, which was adjusted < 0.05)
  mutate(p.adj=round(p.adj, digits = 4)) %>% 
  mutate(y.position = c(y-1.0, y, y-0.5, y-1.4))

#Plot of Day 5 C. diff CFUs with p values added manually----
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
  scale_y_log10(labels=fancy_scientific, breaks = c(10, 100, 10^3, 10^4, 10^5, 10^6, 10^7, 10^8, 10^9), limits = c(NA, 10^9))+
  theme_classic()+
  theme(plot.title=element_text(hjust=0.5))+
  theme(legend.position = "none") + #Get rid of legend 
  theme(text = element_text(size = 16)) +  # Change font size for entire plot
  stat_pvalue_manual(data = pairwise_wilcox_day5_plot, label = "p.adj", y.position = "y.position") 
save_plot(filename = paste0("results/figures/C.diff_CFU_D5_stats.png"), plot_CFU_D5, base_height = 11, base_width = 8.5, base_aspect_ratio = 2)

#Plot of Day 7 C. diff CFUs with p values added manually----
plot_CFU_D7 <- cfu_data_final %>% 
  filter(day == 7) %>% 
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
  scale_y_log10(labels=fancy_scientific, breaks = c(10, 100, 10^3, 10^4, 10^5, 10^6, 10^7, 10^8, 10^9), limits = c(NA, 10^9))+
  theme_classic()+
  theme(plot.title=element_text(hjust=0.5))+
  theme(legend.position = "none") + #Get rid of legend 
  theme(text = element_text(size = 16)) +  # Change font size for entire plot
  stat_pvalue_manual(data = pairwise_wilcox_day7_plot, label = "p.adj", y.position = "y.position") 
save_plot(filename = paste0("results/figures/C.diff_CFU_D7_stats.png"), plot_CFU_D7, base_height = 11, base_width = 8.5, base_aspect_ratio = 2)

#CFU over time plot with astericks on days where cfu varied significantly across sources of mice using annotate()----
#List significant days after BH adjustment of p-values:
#Annotation labels determined by what days were significant
x_annotation <- sig_cfu_days 
y_position <- 10^8.5
label <- cfu_kruskal_stats_adjust %>% 
  mutate(p.value.adj=round(p.value.adj, digits = 4)) %>% 
  filter(p.value.adj <= 0.05) %>%
  mutate(p.signif = case_when(
    p.value.adj <= 0.0001 ~ "****",
    p.value.adj <= 0.001 ~ "***",
    p.value.adj <= 0.01 ~ "**",
    p.value.adj > 0.05 ~ "NS",
    p.value.adj <= 0.05 ~ "*"
  )) %>% 
  pull(p.signif)

#Plot data for all timepoints with annoated significance symbols for test of CFU differences across sources of mice at each timepoint
mean_summary <- cfu_data_final %>% 
  group_by(vendor, day) %>% 
  summarize(mean_cfu = mean(cfu, na.rm = TRUE))
cfu_stats <- ggplot(NULL) +
  geom_point(cfu_data_final, mapping = aes(x = day, y = cfu, color= vendor, fill = vendor), alpha = .2, size = 1.5, show.legend = FALSE, position = position_dodge(width = 0.6)) +
  geom_line(mean_summary, mapping = aes(x = day, y = mean_cfu, color = vendor), size = 1.5) +
  scale_colour_manual(name=NULL,
                      values=color_scheme,
                      breaks=color_vendors,
                      labels=color_vendors)+
  labs(x = "Days Post-Infection", y = "CFU/g Feces") +
  scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9),
                     limits = c(-0.5, 9.5)) +
  geom_hline(yintercept = 100, linetype=2) +
  geom_text(x = 9, y = 104, color = "black", label = "LOD")+
  scale_y_log10(labels=fancy_scientific, breaks = c(10, 100, 10^3, 10^4, 10^5, 10^6, 10^7, 10^8, 10^9), limits = c(NA, 10^9))+
  theme(text = element_text(size = 16))+  # Change font size for entire plot
  annotate("text", y = y_position, x = x_annotation, label = label, size =8)+
  theme_classic() +
  theme(legend.position = c(.92, .8), #(0,0) bottom left (1,1) top right
        legend.key= element_rect(colour = "transparent", fill = "transparent"))
save_plot(filename = paste0("results/figures/cfu_over_time.png"), plot = cfu_stats, base_aspect_ratio = 2)

