source("code/functions.R")

#Graph weight data as percent baseline weight (see Theriot et al. 2011) or weight change (g) from baseline weight.----
#Percent baseline weight for each mouse will be calculated 
#based on the weight recorded on D-1 of the experiment

#Create data frame that has just the baseline_weight data for D-1 of each mouse
baseline_weight <- metadata %>% select(mouse_id, weight, day) %>% 
  filter(day == -1) %>% 
  mutate(baseline_weight = weight) %>% 
  select(mouse_id, baseline_weight)

#Join baseline data frame to main metadata
baseline_weight_data <- inner_join(metadata, baseline_weight, by = "mouse_id") %>% 
  select(experiment, id, mouse_id, vendor, day, weight, baseline_weight) %>% 
  filter(!is.na(weight)) %>% #534 observations that are not NAs
  filter(!day == 10) #Removes 22 observations from Day 10 experiment 2. Drop this timepoint because the data was not collected from any mice from experiment.

#Calculate percent baseline weight data for each mouse based on the D-1 weight.
weight_data <- baseline_weight_data %>%  
  group_by(mouse_id, day) %>% 
  mutate(percent_baseline_weight = 100 + ((weight-baseline_weight)/(baseline_weight))*100) %>% 
  ungroup() %>% 
  group_by(mouse_id) %>% #Group by just mouse_id to figure out the lowest percent baseline weight for each mouse
  mutate(lowest_percent_baseline_weight = min(percent_baseline_weight))  #Create a column to display the lowest percent baseline weight for each mouse

#Calculate weight change (g) for each mouse based on the D-1 weight.
weight_data <- weight_data %>%  
  group_by(mouse_id, day) %>% 
  mutate(weight_change = weight-baseline_weight) %>% 
  ungroup() 

#Function to summarize percent_baseline_weight data (calculate the mean for each group) and plot the data
summarize_plot <- function(df){
  mean_summary <- df %>% 
    group_by(vendor, day) %>% 
    summarize(mean_percent_weight = mean(percent_baseline_weight, na.rm = TRUE))
  ggplot(NULL) +
    geom_point(df, mapping = aes(x = day, y = percent_baseline_weight, color= vendor, fill = vendor), alpha = .2, size = .5, show.legend = FALSE, position = position_dodge(width = 0.6)) +
    geom_line(mean_summary, mapping = aes(x = day, y = mean_percent_weight, color = vendor), size = 1) +
    scale_colour_manual(name=NULL,
                        values=color_scheme,
                        breaks=color_vendors,
                        labels=color_vendors)+
    labs(x = "Days Post-Infection", y = "% Baseline Weight") +
    scale_x_continuous(breaks = c(-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9),
                       limits = c(-1.5, 9.5)) +
    theme_classic()
}

#Combined percent baseline weight plot for the 2 experiments----
combined_exp_percent_weight <- summarize_plot(weight_data)

#Percent baseline weight plot for the 1st experiment----
exp1_percent_weight <- summarize_plot(weight_data %>% filter(experiment == 1))

#Percent baseline weight plot for the 2nd experiment----
exp2_percent_weight <- summarize_plot(weight_data %>% filter(experiment == 2))

plot_grid(combined_exp_percent_weight, exp1_percent_weight, exp2_percent_weight, labels = c("Combined Experiments", "Experiment 1", "Experiment 2"), ncol = 1, label_x = 0, label_y = 1)+
  ggsave("exploratory/notebook/percent_weight_changes.pdf", width = 8.5, height = 11)

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
kruskal_wallis_weight <- weight_data %>% 
  filter(day %in% c('-1', 0, 1, 2, 3, 4, 5, 6, 7, 8, 9)) %>%  #only test days that we have weight data for
  group_by(day) %>% 
  do(tidy(kruskal.test(weight_change~factor(vendor), data=.))) %>% ungroup() %>% 
  mutate(p.value.adj=p.adjust(p.value, method="BH")) %>% 
  arrange(p.value.adj) 
#Timepoints where weight_change is significantly different across the sources of mice
sig_C.diff_weight_timepoints <- kruskal_wallis_cfu %>% 
  filter(p.value.adj <= 0.05) %>% 
  pull(day) 
#Days 5, 6, 7, 4, and 3 are when there are significant differences in weight_change across the different sources of mice (listed in order of increasing adjusted P values)

#For significant timepoints, do pairwise.wilcox.test to determine which sources of mice are significantly different from each other regarding the weight_change.
pairwise.wilcox_groups <- function(timepoint){
  weight_by_day <- weight_data %>% 
    filter(day == timepoint)
  tidy(pairwise.wilcox.test(g = weight_by_day$vendor, x = weight_by_day$weight_change, p.adjust.method = "BH"))
}

# Do pairwise.wilcox tests with BH correction for all significant timepoints----
for(d in sig_C.diff_weight_timepoints){
  name <- paste("pairwise_wilcox_day", d, sep = "") #Way to name the data frames based on the date of interest
  assign(name, pairwise.wilcox_groups(d))
}
pairwise_wilcox_day5 #5 significant: Jackson vs Schloss, Jackson vs Young, Charles River vs Jackson, Taconic vs Jackson, Envigo vs Jackson
pairwise_wilcox_day5 %>% filter(p.value <= 0.05) 
pairwise_wilcox_day6 
pairwise_wilcox_day6 %>% filter(p.value <= 0.05) #0 pairwise comparisons with p < 0.05.
pairwise_wilcox_day7 
pairwise_wilcox_day7 %>% filter(p.value <= 0.05) #4 significant: Jackson vs Schloss, Charles River vs Schloss, Taconic vs Schloss, Jackson vs Young
pairwise_wilcox_day4
pairwise_wilcox_day4 %>% filter(p.value <= 0.05) #0 pairwise comparisons with p < 0.05.
pairwise_wilcox_day3 
pairwise_wilcox_day3 %>% filter(p.value <= 0.05) #0 pairwise comparisons with p < 0.05.

# Boxplots of weight_change data at timepoints where there were significant differences in CFU levels across the different sources of mice:
#Function to plot weight_change data across sources of mice at a specific timepoint----
#Arguments:
# timepoint = timepoint to be analyzed
plot_weight_timepoint <- function(timepoint){
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
    theme_classic()+
    theme(plot.title=element_text(hjust=0.5))+
    theme(legend.position = "none") + #Get rid of legend title & move legend position
    theme(text = element_text(size = 16))  #Remove legend
  save_plot(filename = paste0("results/figures/weight_D", timepoint,".png"), plot_weight_DX, base_height = 11, base_width = 8.5, base_aspect_ratio = 2)
}
#Plot all the timepoints where weight_changes were significantly different across sources of mice
for(d in sig_C.diff_weight_timepoints){
  plot_weight_timepoint(d)
}

