source("code/functions.R")

#Graph weight data as weight change (g) from baseline weight.----

#Range of N mice per day
weight_data %>% group_by(day) %>% count() %>% arrange(n) 

#Plot weight data as Weight change (g)----

#Function to summarize weight change data (calculate the mean for each group) and plot the data
summarize_plot <- function(df){
  mean_summary <- df %>% 
    group_by(vendor, day) %>% 
    summarize(mean_weight_change = mean(weight_change, na.rm = TRUE))
  ggplot(NULL) +
    geom_point(df, mapping = aes(x = day, y = weight_change, color= vendor, fill = vendor), alpha = .2, size = .5, show.legend = FALSE, position = position_dodge(width = 0.6)) +
    geom_line(mean_summary, mapping = aes(x = day, y = mean_weight_change, color = vendor), size = 1) +
    scale_colour_manual(name=NULL,
                        values=color_scheme,
                        breaks=color_vendors,
                        labels=color_vendors)+
    labs(x = "Days Post-Infection", y = "Weight Change (g)") +
    scale_x_continuous(breaks = c(-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9),
                       limits = c(-1.5, 9.5)) +
    theme_classic()
}

#Combined weight change plot for the 2 experiments----
combined_exp_weight <- summarize_plot(weight_data)

#Percent weight change plot for the 1st experiment----
exp1_weight <- summarize_plot(weight_data %>% filter(experiment == 1))

#Percent weight change plot for the 2nd experiment----
exp2_weight <- summarize_plot(weight_data %>% filter(experiment == 2))

plot_grid(combined_exp_weight, exp1_weight, exp2_weight, labels = c("Combined Experiments", "Experiment 1", "Experiment 2"), ncol = 1, label_x = 0, label_y = 1)+
  ggsave("exploratory/notebook/weight_changes.pdf", width = 8.5, height = 11)

#Statistical analysis----
#Kruskal_wallis test for differences across groups at different timepoints with Benjamini-Hochburg correction----
set.seed(19881117) #Same seed used for mothur analysis
weight_kruskal_stats <- weight_data %>% 
  filter(day %in% c('-1', 0, 1, 2, 3, 4, 5, 6, 7, 8, 9)) %>%  #only test days that we have weight data for
  select(day, vendor, weight_change) %>% 
  group_by(day) %>% 
  nest() %>% 
  mutate(model=map(data, ~kruskal.test(x=.x$weight_change, g=as.factor(.x$vendor)) %>% tidy())) %>% 
  mutate(mean = map(data, get_weight_mean_vendor)) %>% 
  unnest(c(model, mean)) %>% 
  ungroup() #Ungroup before adjusting p values

#Adjust p-values for testing multiple days and write results to table
weight_kruskal_stats_adjust <- weight_kruskal_stats %>% 
  select(day, statistic, p.value, parameter, method, Schloss, Young, Jackson, `Charles River`, Taconic, Envigo) %>% 
  mutate(p.value.adj=p.adjust(p.value, method="BH")) %>% 
  arrange(p.value.adj) %>% 
  write_tsv("data/process/weight_stats_all_days.tsv")
#Also write to supplemental table excel file
weight_kruskal_stats_adjust %>% write_xlsx("submission/table_S2_weight_kruskal-wallis.xlsx", format_headers = FALSE)

#Timepoints where weight_change is significantly different across the sources of mice
sig_weight_days <- weight_kruskal_stats_adjust %>% 
  filter(p.value.adj <= 0.05) %>% 
  pull(day) 
#Days 8, 7, 9, 5, 1, 2, 3, and 6 are when there are significant differences in weight_change across the different sources of mice (listed in order of increasing adjusted P values)

#Perform pairwise Wilcoxan rank sum tests for days that were significant by Kruskal-Wallis rank sum test
weight_stats_pairwise <- weight_kruskal_stats %>% 
  filter(day %in% sig_weight_days) %>% #only perform pairwise tests for days that were significant 
  group_by(day) %>% 
  mutate(model=map(data, ~pairwise.wilcox.test(x=.x$weight_change, g=as.factor(.x$vendor), p.adjust.method="BH") %>% 
                     tidy() %>% 
                     mutate(compare=paste(group1, group2, sep="-")) %>% 
                     select(-group1, -group2) %>% 
                     pivot_wider(names_from=compare, values_from=p.value)
  )
  ) %>% 
  unnest(model) %>% 
  select(-data, -parameter, -statistic) %>% 
  write_tsv("data/process/weight_stats_sig_days.tsv")

#Format pairwise stats to use with ggpubr package
plot_format_stats <- weight_stats_pairwise %>%
  #Remove all columns except pairwise comparisons and day
  select(-p.value, -method,-Schloss, -Young, -Jackson, -`Charles River`, -Taconic, -Envigo) %>% 
  group_split() %>% #Keeps a attr(,"ptype") to track prototype of the splits
  lapply(tidy_pairwise) %>% 
  bind_rows()

plot_format_stats %>% filter(day == 8 & p.adj <= 0.05) #7 significant
plot_format_stats %>% filter(day == 7 & p.adj <= 0.05) #4 significant
plot_format_stats %>% filter(day == 9 & p.adj <= 0.05) # 6 significant
plot_format_stats %>% filter(day == 5 & p.adj <= 0.05) #5 significant
plot_format_stats %>% filter(day == 1 & p.adj <= 0.05) #0 significant
plot_format_stats %>% filter(day == 2 & p.adj <= 0.05) #0 significant
plot_format_stats %>% filter(day == 3 & p.adj <= 0.05) #0 significant
plot_format_stats %>% filter(day == 6 & p.adj <= 0.05) # 0 significant

# Boxplots of weight_change data at timepoints where there were significant differences in CFU levels across the different sources of mice:
#Function to plot weight_change data across sources of mice at a specific timepoint----
#Arguments:
# timepoint = timepoint to be analyzed
# stats = data frame of stat values to be added to the plot with ggpubr stat_pvalue_manual
plot_weight_timepoint <- function(timepoint, stats){
  plot_weight_DX <- weight_data %>% 
    filter(day == timepoint) %>% 
    ggplot(aes(x= vendor, y=weight_change, color=vendor))+
    scale_colour_manual(name=NULL,
                        values=color_scheme,
                        breaks=color_vendors,
                        labels=color_vendors)+
    geom_text(x = 11, y = 104, color = "black", label = "LOD")+
    geom_boxplot(outlier.shape = NA, size = 1.2)+
    geom_jitter(shape=19, size=2, alpha=0.6, position=position_jitterdodge(dodge.width=0.7, jitter.width=0.2)) +
    labs(title=NULL, 
         x=NULL,
         y="Weight Change (g)")+
    ylim(-5.5, 3.5)+
    stat_pvalue_manual(data = stats, label = "p.adj", y.position = "y.position")+
    theme_classic()+
    theme(plot.title=element_text(hjust=0.5))+
    theme(legend.position = "none") + #Get rid of legend title & move legend position
    theme(text = element_text(size = 16)) #Remove legend
save_plot(filename = paste0("results/figures/weight_D", timepoint,".png"), plot_weight_DX, base_height = 11, base_width = 8.5, base_aspect_ratio = 2)
}


# Format day 8 stats dataframe
d8_stats <- plot_format_stats %>% 
  filter(day == 8 & p.adj <= 0.05) %>% 
  mutate(p.adj=round(p.adj, digits = 4)) %>% 
  mutate(y.position = c(1.4, 2.6, 3.4, 2, 3, 2.3, 1.8))
#Plot of day 8----
plot_weight_timepoint(8, d8_stats)
# Format day 7 stats dataframe
d7_stats <- plot_format_stats %>% 
  filter(day == 7 & p.adj <= 0.05) %>% 
  mutate(p.adj=round(p.adj, digits = 4)) %>% 
  mutate(y.position = c(2.8, 3.2, 3.5, 2.4))
#Plot of day 7----
plot_weight_timepoint(7, d7_stats)
# Format day 9 stats dataframe
d9_stats <- plot_format_stats %>% 
  filter(day == 9 & p.adj <= 0.05) %>% 
  mutate(p.adj=round(p.adj, digits = 4)) %>% 
  mutate(y.position = c(1.4, 2.6, 3.4, 2.2, 3, 1.8))
# Plot of day 9---- 
plot_weight_timepoint(9, d9_stats)
# Format day 5 stats dataframe
d5_stats <- plot_format_stats %>% 
  filter(day == 5 & p.adj <= 0.05) %>% 
  mutate(p.adj=round(p.adj, digits = 4)) %>% 
  mutate(y.position = c(2.9, 2.2, 2.6, 3.2, 3.5))
# Plot of day 5----
plot_weight_timepoint(5, d5_stats)
#Day 1,2,3, and 6 have no significant pairwise comparisons, so there's no need to format the stats into a data frame.
#Create an empty place holder data frame to use as 2nd argument for plot_weight_timepoints function for days 1,2,3 and 6
stats_na <- plot_format_stats %>% 
  filter(day == 1 & p.adj <= 0.05) %>% 
  mutate(p.adj=round(p.adj, digits = 4)) %>% 
  mutate(y.position = NA)
# Plot of day 1----
plot_weight_timepoint(1, stats_na)
#Plot of day 2---- 
plot_weight_timepoint(2, stats_na)
#Plot of day 3
plot_weight_timepoint(3, stats_na)
#Plot of day 6
plot_weight_timepoint(6, stats_na)

#Weight over time plot with astericks on days where weight varied significantly across sources of mice using annotate()----
#List significant days after BH adjustment of p-values:
#Annotation labels determined by what days were significant
x_annotation <- sig_weight_days 
y_position <- 3.4
label <- weight_kruskal_stats_adjust %>% 
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

mean_summary <- weight_data %>% 
    group_by(vendor, day) %>% 
    summarize(mean_weight_change = mean(weight_change, na.rm = TRUE))
weight_stats <-  ggplot(NULL) +
    geom_point(weight_data, mapping = aes(x = day, y = weight_change, color= vendor, fill = vendor), alpha = .2, size = 1.5, show.legend = FALSE, position = position_dodge(width = 0.6)) +
    geom_line(mean_summary, mapping = aes(x = day, y = mean_weight_change, color = vendor), size = 1.5) +
    scale_colour_manual(name=NULL,
                        values=color_scheme,
                        breaks=color_vendors,
                        labels=color_vendors)+
    labs(x = "Days Post-Infection", y = "Weight Change (g)") +
    scale_x_continuous(breaks = c(-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9),
                       limits = c(-1.5, 9.5)) +
    theme(text = element_text(size = 16))+  # Change font size for entire plot
    annotate("text", y = y_position, x = x_annotation, label = label, size =8)+
    theme_classic()
save_plot(filename = paste0("results/figures/weight_over_time.png"), plot = weight_stats, base_aspect_ratio = 2)


