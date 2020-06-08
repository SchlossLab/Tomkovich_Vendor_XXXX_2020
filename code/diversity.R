source("code/functions.R")

diversity_data <- read_tsv("data/process/vendors.diversity.summary") %>%
  filter(method == "ave") %>% 
  select(group, sobs, shannon, invsimpson, coverage) %>% 
  inner_join(metadata, by = c("group" = "id")) %>% #Match only the samples we have sequence data for
  mutate(day=factor(day, levels=c("-1", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9"))) # Treat day as a factor

#Boxplots of diversity metrics (shannon and richness) for all sources of mice at all timepoints:
shannon_all <- diversity_data %>% 
  ggplot(aes(x=vendor, y = shannon, colour= vendor))+
  scale_colour_manual(name=NULL,
                      values=color_scheme,
                      breaks=color_vendors,
                      labels=color_vendors)+
  geom_boxplot(outlier.shape = NA, size = 1.2)+
  geom_jitter(shape=19, size=2, alpha=0.6, position=position_jitterdodge(dodge.width=0.7, jitter.width=0.1), show.legend = FALSE) +
  labs(title=NULL, 
       x="Source",
       y="Shannon Diversity Index")+
  theme_classic()+
  theme(legend.position = "none",
        text = element_text(size = 16))+ # Change font size for entire plot
  facet_wrap(vars(day))

richness_all <- diversity_data %>% 
  ggplot(aes(x=vendor, y = sobs, colour= vendor))+
  scale_colour_manual(name=NULL,
                      values=color_scheme,
                      breaks=color_vendors,
                      labels=color_vendors)+
  geom_boxplot(outlier.shape = NA, size = 1.2)+
  geom_jitter(shape=19, size=2, alpha=0.6, position=position_jitterdodge(dodge.width=0.7, jitter.width=0.1), show.legend = FALSE) +
  labs(title=NULL, 
       x="Source",
       y="Number of Observed OTUs")+
  theme_classic()+
  theme(legend.position = "none",
        text = element_text(size = 16))+ # Change font size for entire plot
  facet_wrap(vars(day))

#Statistical analysis
set.seed(19881117) #Same seed used for mothur analysis

#Test if Shannon and richness are normally distributed with Shapiro-Wilk normality test----
#Shapiro-Wilk test of shannon data
shannon_by_day <- diversity_data %>% 
  group_by(day) %>% 
  nest() %>% 
  mutate(Shapiro = map(data, ~shapiro.test(.x$shannon)))

shannon_by_day_glance <- shannon_by_day %>%
  mutate(glance_shapiro = Shapiro %>% map(glance)) %>%
  unnest(glance_shapiro)
#p < 0.05 means the data are not normally distributed
# Days 2, 0 are not normally distrubuted. Shannon data is normally distributed for the rest of the timepoints

#Shapiro-Wilk test of richness data
richness_by_day <- diversity_data %>% 
  group_by(day) %>% 
  nest() %>% 
  mutate(Shapiro = map(data, ~shapiro.test(.x$sobs)))

richness_by_day_glance <- richness_by_day %>%
  mutate(glance_shapiro = Shapiro %>% map(glance)) %>%
  unnest(glance_shapiro)
#p < 0.05 means the data are not normally distributed
#For all timepoints, richness data is not normally distributed

#Since the majority of diversity data is not normally distributed will test for differences across sources of mice at specific timepoints with the Kruskal-Wallis test 
#Function for Kruskal_wallis test for differences in Shannon index across groups with Benjamini-Hochburg correction----
shannon_stats <- diversity_data %>% 
    select(vendor, shannon, day) %>% 
    group_by(day) %>% 
    nest() %>% 
    mutate(model=map(data, ~kruskal.test(x=.x$shannon, g=as.factor(.x$vendor)) %>% tidy())) %>% 
    mutate(mean = map(data, get_shannon_mean_vendor)) %>% 
    unnest(c(model, mean)) %>% 
    ungroup() 

#Adjust p-values for testing multiple days and write results to table
shannon_kruskal_stats_adjust <- shannon_stats %>% 
  select(day, statistic, p.value, parameter, method, Schloss, Young, Jackson, `Charles River`, Taconic, Envigo) %>% 
  mutate(p.value.adj=p.adjust(p.value, method="BH")) %>% 
  arrange(p.value.adj) %>% 
  write_tsv("data/process/shannon_stats_all_days.tsv")

#List significant days after BH adjustment of p-values:
sig_shannon_days <- shannon_kruskal_stats_adjust %>% 
  filter(p.value.adj <= 0.05) %>% 
  pull(day)
# There were significant differences in Shannon index across vendors on days 1, 5, 9, 3, 0, 6 and 4 

#Perform pairwise Wilcoxan rank sum tests for days that were significant by Kruskal-Wallis rank sum test
shannon_stats_pairwise <- shannon_stats %>% 
  filter(day %in% sig_shannon_days) %>% #only perform pairwise tests for days that were significant 
  group_by(day) %>% 
  mutate(model=map(data, ~pairwise.wilcox.test(x=.x$shannon, g=as.factor(.x$vendor), p.adjust.method="BH") %>% 
                     tidy() %>% 
                     mutate(compare=paste(group1, group2, sep="-")) %>% 
                     select(-group1, -group2) %>% 
                     pivot_wider(names_from=compare, values_from=p.value)
  )
  ) %>% 
  unnest(model) %>% 
  select(-data, -parameter, -statistic, -p.value, -method, -Schloss, -Young, -Jackson, -`Charles River`, -Taconic, -Envigo) %>% 
#Combine with shannon_kruskal_stats_adjust so that adjusted p-values are on the same table
  inner_join(shannon_kruskal_stats_adjust, by = c("day")) %>% 
  select(-p.value, -parameter, -statistic, -Schloss, -Young, -Jackson, -`Charles River`, -Taconic, -Envigo) %>% 
  write_tsv("data/process/shannon_stats_sig_days.tsv")

#Format shannon pairwise stats to use with ggpubr package
plot_format_shannon <- shannon_stats_pairwise %>%
  #Remove all columns except pairwise comparisons and day
  select(-method, -p.value.adj) %>% 
  group_split() %>% #Keeps a attr(,"ptype") to track prototype of the splits
  lapply(tidy_pairwise) %>% 
  bind_rows()

#Function for Kruskal_wallis test for differences in sobs (richness) across groups with Benjamini-Hochburg correction----
richness_stats <- diversity_data %>% 
  select(vendor, sobs, day) %>% 
  group_by(day) %>% 
  nest() %>% 
  mutate(model=map(data, ~kruskal.test(x=.x$sobs, g=as.factor(.x$vendor)) %>% tidy())) %>% 
  mutate(mean = map(data, get_sobs_mean_vendor)) %>% 
  unnest(c(model, mean)) %>% 
  ungroup() 

#Adjust p-values for testing multiple days and write results to table
richness_kruskal_stats_adjust <- richness_stats %>% 
  select(day, statistic, p.value, parameter, method, Schloss, Young, Jackson, `Charles River`, Taconic, Envigo) %>% 
  mutate(p.value.adj=p.adjust(p.value, method="BH")) %>% 
  arrange(p.value.adj) %>% 
  write_tsv("data/process/richness_stats_all_days.tsv")

#Combine Kruskal-Wallis test results for shannon and richness diversity metrics
shannon_results <- shannon_kruskal_stats_adjust %>% 
  mutate(diversity_metric = "Shannon Diversity")
richness_results <- richness_kruskal_stats_adjust %>% 
  mutate(diversity_metric = "Richness")
diversity_kruskal_stats_adjust <- rbind(shannon_results, richness_results) %>% 
  write_xlsx("submission/table_S3_diversity_kruskal-wallis.xlsx", format_headers = FALSE)

#List significant days after BH adjustment of p-values:
sig_richness_days <- richness_kruskal_stats_adjust %>% 
  filter(p.value.adj <= 0.05) %>% 
  pull(day)
# There were significant differences in richness (sobs) across vendors on days 5, 9, 6, 1, 3, 7, 4, 0, and -1 

#Perform pairwise Wilcoxan rank sum tests for days that were significant by Kruskal-Wallis rank sum test
richness_stats_pairwise <- richness_stats %>% 
  filter(day %in% sig_richness_days) %>% #only perform pairwise tests for days that were significant 
  group_by(day) %>% 
  mutate(model=map(data, ~pairwise.wilcox.test(x=.x$sobs, g=as.factor(.x$vendor), p.adjust.method="BH") %>% 
                     tidy() %>% 
                     mutate(compare=paste(group1, group2, sep="-")) %>% 
                     select(-group1, -group2) %>% 
                     pivot_wider(names_from=compare, values_from=p.value)
  )
  ) %>% 
  unnest(model) %>% 
  select(-data, -parameter, -statistic, -p.value, -method, -Schloss, -Young, -Jackson, -`Charles River`, -Taconic, -Envigo) %>% 
#Combine with richness_kruskal_stats_adjust so that adjusted p-values are on the same table
  inner_join(richness_kruskal_stats_adjust, by = c("day")) %>% 
  select(-p.value, -parameter, -statistic, -Schloss, -Young, -Jackson, -`Charles River`, -Taconic, -Envigo) %>% 
  write_tsv("data/process/richness_stats_sig_days.tsv")

#Combine pairwise wilcoxon test results for shannon and richness diversity metrics
shannon_pairwise_results <- shannon_stats_pairwise %>% 
  mutate(diversity_metric = "Shannon Diversity")
richness_pairwise_results <- richness_stats_pairwise %>% 
  mutate(diversity_metric = "Richness")
diversity_pairwise_results <- rbind(shannon_pairwise_results, richness_pairwise_results) %>% 
  write_xlsx("submission/table_S4_diversity_pairwise.xlsx", format_headers = FALSE)

#Format richness pairwise stats to use with ggpubr package
plot_format_richness <- richness_stats_pairwise %>%
  #Remove all columns except pairwise comparisons and day
  select(-method,-p.value.adj) %>% 
  group_split() %>% #Keeps a attr(,"ptype") to track prototype of the splits
  lapply(tidy_pairwise) %>% 
  bind_rows()

#Function to make boxplots of shannon for all sources of mice on a specific experimental day----
#Arguments: 
#  timepoint = day of the experiment
shannon_dx_plot <- function(timepoint){
  diversity_data %>% 
    filter(day == timepoint) %>% 
    ggplot(aes(x=vendor, y =shannon, colour= vendor))+
    scale_colour_manual(name=NULL,
                        values=color_scheme,
                        breaks=color_vendors,
                        labels=color_vendors)+
    scale_x_discrete(guide = guide_axis(n.dodge = 2))+
    geom_boxplot(outlier.shape = NA, size = 1.2)+
    geom_jitter(shape=19, size=2, alpha=0.6, position=position_jitterdodge(dodge.width=0.7, jitter.width=0.1), show.legend = FALSE) +
    labs(title=NULL, 
         x="Source",
         y="Shannon Diversity Index")+
    ylim(0, 5.5)+
    theme_classic()+
    theme(legend.position = "none",
          text = element_text(size = 16)) # Change font size for entire plot
}

#Boxplots of shannon for all sources of mice on day -1, the baseline microbiota community for each mouse----
shannon_dn1 <- shannon_dx_plot(-1) +
  ggtitle("Baseline")+ #Title plot
  theme(plot.title = element_text(hjust = 0.5)) #Center plot titile
#There was no overall significant difference in Shannon index at this timepoint across sources of mice.  
save_plot("results/figures/shannon_dn1.png", shannon_dn1) #Use save_plot instead of ggsave because it works better with cowplot

#Add p.value manually to timepoints of interest with ggpubr stat_pvalue_manual() function: ----
#Data frames of day 0 p.values to add manually
pairwise_shannon_day0_plot <- plot_format_shannon %>% 
  filter(day == 0) %>%  
  filter(p.adj <= 0.05) %>% #Only show comparisons that were significant. p.value, which was adjusted < 0.05)
  mutate(p.adj=round(p.adj, digits = 4)) %>% 
  mutate(y.position = c(3.5, 4, 3, 2.5))
#4 pairwise comparisons were significant at day 0. 
#Day 0 Plot with stats for pairwise comparisons:
shannon_d0 <- shannon_dx_plot(0) +
  ggtitle("Clindamycin")+ #Title plot
  theme(plot.title = element_text(hjust = 0.5)) +#Center plot titile
  stat_pvalue_manual(data = pairwise_shannon_day0_plot, label = "p.adj", y.position = "y.position") 
save_plot("results/figures/shannon_d0.png", shannon_d0)

#Data frames of day 1 p.values to add manually
pairwise_shannon_day1_plot <- plot_format_shannon %>% 
  filter(day == 1) %>%  
  filter(p.adj <= 0.05) %>% #Only show comparisons that were significant. p.value, which was adjusted < 0.05)
  mutate(p.adj=round(p.adj, digits = 4)) %>% 
  mutate(y.position = c(4.7, 5.3, 3.8, 5, 3.1, 4.2, 4.5, 3.55))
#8 pairwise comparisons were significant at day 0. 
shannon_d1 <- shannon_dx_plot(1) +
  ggtitle("Post-infection")+ #Title plot
  theme(plot.title = element_text(hjust = 0.5)) +#Center plot titile
  stat_pvalue_manual(data = pairwise_shannon_day1_plot, label = "p.adj", y.position = "y.position") 
save_plot("results/figures/shannon_d1.png", shannon_d1)

#Function to make boxplots of sobs (richness) for all sources of mice on a specific experimental day----
#Arguments: 
#  timepoint = day of the experiment
sobs_dx_plot <- function(timepoint){
  diversity_data %>% 
  filter(day == timepoint) %>% 
  ggplot(aes(x=vendor, y =sobs, colour= vendor))+
  scale_colour_manual(name=NULL,
                      values=color_scheme,
                      breaks=color_vendors,
                      labels=color_vendors)+
  geom_boxplot(outlier.shape = NA, size = 1.2)+
  scale_x_discrete(guide = guide_axis(n.dodge = 2))+
  geom_jitter(shape=19, size=2, alpha=0.6, position=position_jitterdodge(dodge.width=0.7, jitter.width=0.1), show.legend = FALSE) +
  labs(title=NULL, 
       x="Source",
       y="Number of Observed OTUs")+
  ylim(0, 320)+
  theme_classic()+
  theme(legend.position = "none",
        text = element_text(size = 16)) # Change font size for entire plot
}
  
#Boxplots of sobs (richness) for all sources of mice on day -1, the baseline microbiota community for each mouse----
#Add p.value manually to timepoints of interest with ggpubr stat_pvalue_manual() function: ----
#Data frame of day -1 p.values to add manually
pairwise_sobs_dayn1_plot <- plot_format_richness %>% 
  filter(day == -1) %>%  
  filter(p.adj <= 0.05)  #Only show comparisons that were significant. p.value, which was adjusted < 0.05)
#No pairwise comparisons were significant at day -1.  
# Plot of richness across sources of mice on day -1
sobs_dn1 <- sobs_dx_plot(-1)+
  ggtitle("Baseline")+ #Title plot
  theme(plot.title = element_text(hjust = 0.5)) #Center plot titile
save_plot("results/figures/richness_dn1.png", sobs_dn1) #Use save_plot instead of ggsave because it works better with cowplot

#Data frame of day 0 p.values to add manually
pairwise_sobs_day0_plot <- plot_format_richness %>% 
  filter(day == 0) %>%  
  filter(p.adj <= 0.05) %>%  #Only show comparisons that were significant. p.value, which was adjusted < 0.05)
  mutate(p.adj=round(p.adj, digits = 4)) %>% 
  mutate(y.position = c(150, 200, 250))
  #3 significant pairwise
# Plot of richness across sources of mice on day 0 with significant pairwise comparison p values
sobs_d0 <- sobs_dx_plot(0) +
  stat_pvalue_manual(data = pairwise_sobs_day0_plot, label = "p.adj", y.position = "y.position") +
  ggtitle("Clindamycin")+ #Title plot
  theme(plot.title = element_text(hjust = 0.5)) #Center plot titile
save_plot("results/figures/richness_d0.png", sobs_d0) #Use save_plot instead of ggsave because it works better with cowplot

#Data frame of day 1 p.values to add manually
pairwise_sobs_day1_plot <- plot_format_richness %>% 
  filter(day == 1) %>%  
  filter(p.adj <= 0.05) %>%  #Only show comparisons that were significant. p.value, which was adjusted < 0.05)
  mutate(p.adj=round(p.adj, digits = 4)) %>% 
  mutate(y.position = c(250, 200, 150, 300))
#4 significant pairwise

# Plot of richness across sources of mice on day 1 with significant pairwise comparison p values
sobs_d1 <- sobs_dx_plot(1) +
  stat_pvalue_manual(data = pairwise_sobs_day1_plot, label = "p.adj", y.position = "y.position") +
  ggtitle("Post-infection")+ #Title plot
  theme(plot.title = element_text(hjust = 0.5)) #Center plot titile
save_plot("results/figures/richness_d1.png", sobs_d1) #Use save_plot instead of ggsave because it works better with cowplot

