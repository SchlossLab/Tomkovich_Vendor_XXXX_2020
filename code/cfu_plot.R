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
  mutate(mean = map(data, get_cfu_mean_vendor)) %>% 
  unnest(c(model, mean)) %>% 
  ungroup() #Ungroup before adjusting p values

#Adjust p-values for testing multiple days and write results to table
cfu_kruskal_stats_adjust <- cfu_kruskal_stats %>% 
  select(day, statistic, p.value, parameter, method, Schloss, Young, Jackson, `Charles River`, Taconic, Envigo) %>% 
  mutate(p.value.adj=p.adjust(p.value, method="BH")) %>% 
  arrange(p.value.adj) %>% 
  write_tsv("data/process/cfu_stats_all_days.tsv")
#Also write to supplemental table excel file
cfu_kruskal_stats_adjust %>% write_xlsx("submission/table_S1_cfu_kruskal-wallis.xlsx", format_headers = FALSE)

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
    geom_text(x = 9.2, y = 102, color = "black", label = "LOD")+
    scale_y_log10(labels=fancy_scientific, breaks = c(10, 100, 10^3, 10^4, 10^5, 10^6, 10^7, 10^8, 10^9), limits = c(40, 350000000))+
    theme_classic()
}

#CFU plot for the 1st experiment----
exp1_cfu <- summarize_plot(cfu_data_final %>% filter(experiment == 1))+
  ggtitle("1st experiment")+ #Title plot
  theme(plot.title = element_text(hjust = 0.5)) #Center plot titile
save_plot(filename = paste0("results/figures/exp1_cfu_time.png"), plot = exp1_cfu, base_aspect_ratio = 2)

#CFU plot for the 2nd experiment----
exp2_cfu <- summarize_plot(cfu_data_final %>% filter(experiment == 2))+
  ggtitle("2nd experiment")+ #Title plot
  theme(plot.title = element_text(hjust = 0.5)) #Center plot titile
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
    p.value.adj <= 0.0001 ~ "****",
    p.value.adj <= 0.001 ~ "***",
    p.value.adj <= 0.01 ~ "**",
    p.value.adj > 0.05 ~ "NS",
    p.value.adj <= 0.05 ~ "*"
  )) %>% 
  pull(p.signif)

#Plot data for all timepoints with annotated significance symbols for test of CFU differences across sources of mice at each timepoint
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

#Plots of day 7 colonization status input data used to build the logistic regression classification models

#For day -1 mice with sequence data
dn1_classification_ids <- read_csv("data/process/classification_input_ids_dn1.csv") %>% pull(day_n1)
d7_title <- c(expression(paste("Day 7 ", italic("C. difficile"), " status for day -1 input data"))) #Expression variable for the title so that bacteria name will be in italics
d7_status_dn1_input <- metadata %>% #Use metadata because day-1 did not have any CFU data for it and is filtered out of cfu_data_final
  filter(id %in% dn1_classification_ids) %>% 
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
save_plot(filename = paste0("results/figures/d7_status_dn1_input.png"), plot = d7_status_dn1_input, base_aspect_ratio = 2)

#For day 0 mice with sequence data
d0_classification_ids <- read_csv("data/process/classification_input_ids_d0.csv") %>% pull(day_0)
d7_title <- c(expression(paste("Day 7 ", italic("C. difficile"), " status for day 0 input data"))) #Expression variable for the title so that bacteria name will be in italics
d7_status_d0_input <- metadata %>% #Use metadata because day-1 did not have any CFU data for it and is filtered out of cfu_data_final
  filter(id %in% d0_classification_ids) %>% 
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
save_plot(filename = paste0("results/figures/d7_status_d0_input.png"), plot = d7_status_d0_input, base_aspect_ratio = 2)

#For day 1 mice with sequence data
d1_classification_ids <- read_csv("data/process/classification_input_ids_d1.csv") %>% pull(day_1)
d7_title <- c(expression(paste("Day 7 ", italic("C. difficile"), " status for day 1 input data"))) #Expression variable for the title so that bacteria name will be in italics
d7_status_d1_input <- metadata %>% #Use metadata because day-1 did not have any CFU data for it and is filtered out of cfu_data_final
  filter(id %in% d1_classification_ids) %>% 
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
save_plot(filename = paste0("results/figures/d7_status_d1_input.png"), plot = d7_status_d1_input, base_aspect_ratio = 2)

