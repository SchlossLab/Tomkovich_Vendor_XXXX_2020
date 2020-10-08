source("code/functions.R")

#Graph weight data as weight change (g) from baseline weight.----
weight_data <- metadata %>%
  filter(!is.na(weight)) #512 observations that are not NAs

#Range of N mice per day
weight_data %>% group_by(day) %>% count() %>% arrange(n) 

#Plot weight data as Weight change (g)----

#Function to summarize weight change data (calculate the median for each group) and plot the data
summarize_plot <- function(df){
  median_summary <- df %>% 
    group_by(vendor, day) %>% 
    summarize(median_weight_change = median(weight_change, na.rm = TRUE))
  ggplot(NULL) +
    geom_point(df, mapping = aes(x = day, y = weight_change, color= vendor, fill = vendor, shape = experiment), size = .5, show.legend = FALSE, position = position_dodge(width = 0.6)) +
    geom_line(median_summary, mapping = aes(x = day, y = median_weight_change, color = vendor), size = 1) +
    scale_colour_manual(name=NULL,
                        values=color_scheme,
                        breaks=color_vendors,
                        labels=color_vendors)+
    scale_shape_manual(name=NULL,
                       values=shape_scheme,
                       breaks=shape_experiment,
                       labels=shape_experiment) +
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

#Statistical analysis----
#Kruskal_wallis test for differences across groups at different timepoints with Benjamini-Hochburg correction----
set.seed(19881117) #Same seed used for mothur analysis
weight_kruskal_stats <- weight_data %>% 
  filter(day %in% c('-1', 0, 1, 2, 3, 4, 5, 6, 7, 8, 9)) %>%  #only test days that we have weight data for
  select(day, vendor, weight_change) %>% 
  group_by(day) %>% 
  nest() %>% 
  mutate(model=map(data, ~kruskal.test(x=.x$weight_change, g=as.factor(.x$vendor)) %>% tidy())) %>% 
  mutate(median = map(data, get_weight_median_vendor)) %>% 
  unnest(c(model, median)) %>% 
  ungroup() #Ungroup before adjusting p values

#Adjust p-values for testing multiple days and write results to table
weight_kruskal_stats_adjust <- weight_kruskal_stats %>% 
  select(day, statistic, p.value, parameter, method, Schloss, Young, Jackson, `Charles River`, Taconic, Envigo) %>% 
  mutate(p.value.adj=p.adjust(p.value, method="BH")) %>% 
  arrange(p.value.adj) %>% 
  write_tsv("data/process/weight_stats_all_days.tsv")

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

#Weight over time plot with astericks on days where weight varied significantly across sources of mice using annotate()----
#List significant days after BH adjustment of p-values:
#Annotation labels determined by what days were significant
x_annotation <- sig_weight_days 
y_position <- 3.4
label <- weight_kruskal_stats_adjust %>% 
  mutate(p.value.adj=round(p.value.adj, digits = 4)) %>% 
  filter(p.value.adj <= 0.05) %>%
  mutate(p.signif = case_when(
    p.value.adj > 0.05 ~ "NS",
    p.value.adj <= 0.05 ~ "*"
  )) %>% 
  pull(p.signif)

median_summary <- weight_data %>% 
    group_by(vendor, day) %>% 
    summarize(median_weight_change = median(weight_change, na.rm = TRUE))
weight_stats <-  ggplot(NULL) +
    geom_point(weight_data, mapping = aes(x = day, y = weight_change, color= vendor, fill = vendor, shape = experiment), size = 1.5, show.legend = FALSE, position = position_dodge(width = 0.6)) +
    geom_line(median_summary, mapping = aes(x = day, y = median_weight_change, color = vendor), size = 1.5) +
    scale_colour_manual(name=NULL,
                        values=color_scheme,
                        breaks=color_vendors,
                        labels=color_vendors)+
    scale_shape_manual(name=NULL,
                        values=shape_scheme,
                        breaks=shape_experiment,
                        labels=shape_experiment) +
    labs(x = "Days Post-Infection", y = "Weight Change (g)") +
    scale_x_continuous(breaks = c(-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9),
                       limits = c(-1.5, 9.5)) +
    theme(text = element_text(size = 16))+  # Change font size for entire plot
    annotate("text", y = y_position, x = x_annotation, label = label, size =8)+
    theme_classic()+
    theme(legend.position = "none", #Get rid of legend 
          text = element_text(size = 18),
          panel.grid.minor.x = element_line(size = 0.4, color = "grey"))  # Add gray lines to clearly separate symbols by days)
save_plot(filename = paste0("results/figures/weight_over_time.png"), plot = weight_stats, base_aspect_ratio = 2)

#Add p.value manually to timepoints of interest with ggpubr stat_pvalue_manual() function: 2 ----
#Figure out y.position via log transformation to match y-axis scale
y <- max(weight_data$weight_change)

pairwise_wilcox_day2_plot <- plot_format_stats %>% 
  filter(day == 2) %>% 
  filter(p.adj <= 0.05)  #Only show comparisons that were significant. p.value, which was adjusted < 0.05)
#No pairwise comparisons were significant on day 2

#Plot of Day 2 weight change data with p values added manually----
median_d2 <- weight_data %>% 
  filter(day == 2) %>% 
  group_by(vendor) %>% 
  mutate(median_weight = median(weight_change, na.rm = TRUE)) %>% #create a column of median_cfu
  ungroup()
plot_weight_D2 <- median_d2 %>% 
  ggplot(aes(x= vendor, y=weight_change, color=vendor))+
  scale_colour_manual(name=NULL,
                      values=color_scheme,
                      breaks=color_vendors,
                      labels=color_vendors)+
  geom_errorbar(aes(ymax = median_weight, ymin = median_weight), color = "gray50", size = 2)+ #Add lines to indicate the median for each group to the plot. Median calculated before y axis transformation
  geom_jitter(aes(shape=experiment, size=2)) +
  scale_shape_manual(name=NULL,
                     values=shape_scheme,
                     breaks=shape_experiment,
                     labels=shape_experiment) +
  labs(title = "Day 2 Post-infection", x = NULL, y = "Weight Change (g)") +
  scale_x_discrete(guide = guide_axis(n.dodge = 2))+  
  theme_classic()+
  theme(plot.title=element_text(hjust=0.5))+
  theme(legend.position = "none") + #Get rid of legend 
  theme(text = element_text(size = 16)) # Change font size for entire plot
save_plot(filename = paste0("results/figures/weight_D2.png"), plot_weight_D2, base_height = 11, base_width = 8.5, base_aspect_ratio = 2)

