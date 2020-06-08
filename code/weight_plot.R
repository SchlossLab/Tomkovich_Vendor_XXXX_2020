source("code/functions.R")

#Graph weight data as weight change (g) from baseline weight.----
weight_data <- metadata %>%
  filter(!is.na(weight)) #512 observations that are not NAs

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


