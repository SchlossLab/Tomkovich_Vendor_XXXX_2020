source("code/functions.R")

read_dist_df <- function(dist_file_name){
  linear_data <- scan(dist_file_name, what="character", quiet=TRUE)[-1]
  
  samples <- str_subset(linear_data, "D") #Pull out sample ids with E (denotes experiment number
  n_samples <- length(samples)
  distance_strings <- str_subset(linear_data, "\\.")
  
  distance_matrix <- matrix(0, nrow=n_samples, ncol=n_samples)
  colnames(distance_matrix) <- samples
  as_tibble(cbind(rows=samples, distance_matrix)) %>%
    gather(columns, distances, -rows) %>%
    filter(rows < columns) %>%
    arrange(columns, rows) %>%
    mutate(distances = as.numeric(distance_strings))
}

#Day -1 Inter and Intra-group variation----
dn1_dist <- read_dist_df("data/mothur/d-1/vendors.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.0.03.subsample.thetayc.0.03.lt.ave.dist") %>% 
  #Relabel rows and columns by separating sample id into mouse and exp labels
  separate(rows, into=c("row_mouse", "row_day", "row_exp"), sep = c(3, -2)) %>% #separate mouse ID from D and experiment number
  separate(columns, into=c("col_mouse", "col_day", "col_exp"), sep = c(3, -2)) %>% #separate mouse ID from D and experiment number
  select(-row_day, -col_day) %>% #get rid of day column, since all samples are from day -1 timepoint
  mutate(experiment = recode(row_exp, E1 = 1, E2 = 2)) %>% #Match experiment column variable names from metadata
  mutate(col_exp = recode(col_exp, E1 = 1, E2 = 2)) %>% #Match experiment column variable names from metadata
  #Make 2 columns to designate row and column vendors
  mutate(vendor = case_when(str_detect(row_mouse, "S") ~ "Schloss",
                                str_detect(row_mouse, "Y") ~ "Young",
                                str_detect(row_mouse, "J") ~ "Jackson",
                                str_detect(row_mouse, "C") ~ "Charles River",
                                str_detect(row_mouse, "T") ~ "Taconic",
                                str_detect(row_mouse, "E") ~ "Envigo")
         ) %>% 
  mutate(col_vendor = case_when(str_detect(col_mouse, "S") ~ "Schloss",
                                str_detect(col_mouse, "Y") ~ "Young",
                                str_detect(col_mouse, "J") ~ "Jackson",
                                str_detect(col_mouse, "C") ~ "Charles River",
                                str_detect(col_mouse, "T") ~ "Taconic",
                                str_detect(col_mouse, "E") ~ "Envigo")
  ) %>% 
  mutate(experiment=factor(experiment, levels=c("1", "2")), # Make sure experiment is treated as a factor
         vendor=factor(vendor, levels=c("Schloss", "Young", "Jackson", "Charles River", "Taconic", "Envigo"))) #Transform vendor into factor, same as what was done for metadata

#Within experiment and same vendor comparisons
within_exp <- dn1_dist %>% 
  filter(vendor == col_vendor) %>% 
  filter(experiment == col_exp) %>% 
  group_by(vendor) %>% 
#  summarize(median_dist = median(distances), n = n()) %>% 
  ggplot(aes(x = vendor, y = distances, color = vendor)) +
  geom_boxplot(outlier.shape = NA, size = 1.2, show.legend = FALSE)+
  geom_jitter(aes(shape = experiment), size=2, alpha=0.6, show.legend = FALSE) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2))+
  ylim(0, 1.25)+
  scale_colour_manual(name=NULL,
                      values=color_scheme,
                      breaks=color_vendors,
                      labels=color_vendors)+
  labs(title = "Within experiment", x = NULL, y = "Theta YC Distance")+
  theme(plot.title = element_text(hjust = 0.5)) +#Center plot title
  theme_classic()
save_plot("results/figures/within_exp_dn1.png", within_exp)

#Between experiment and same vendor comparisons
between_exp <- dn1_dist %>% 
  filter(vendor == col_vendor) %>% 
  filter(experiment != col_exp) %>% 
  group_by(vendor) %>% 
  #  summarize(median_dist = median(distances), n = n()) %>% 
  ggplot(aes(x = vendor, y = distances, color = vendor)) +
  geom_boxplot(outlier.shape = NA, size = 1.2)+
  geom_jitter(aes(shape = experiment), size=2, alpha=0.6, show.legend = FALSE) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2))+
  ylim(0, 1.25)+
  scale_colour_manual(name=NULL,
                      values=color_scheme,
                      breaks=color_vendors,
                      labels=color_vendors)+
  labs(title = "Between experiments", x = NULL, y = "Theta YC Distance")+
  theme(plot.title = element_text(hjust = 0.5))+ #Center plot title
  theme_classic()
save_plot("results/figures/between_exp_dn1.png", between_exp)

#Day -1 versus D7 thetayc distances in individual mice, compared across sources----
#read in theta YC distances for all samples, all timepoints
all_dist <- read_dist_df("data/process/vendors.subsample.thetayc.ave.dist")

#Modify metadata to match row_unique_mouse, select clearance_status_d7 column and join to dn1vd7_dist
mouse_clear_d7_status <- metadata %>% 
  separate(id, into=c("row_mouse", "comparison_day", "row_exp"), sep = c(3, -2)) %>% #separate mouse ID from D and experiment number
  unite("row_unique_mouse", sep = "_", row_mouse, row_exp, remove = TRUE) %>% #make columns to denote unique individual mouse from columns
  distinct(row_unique_mouse, .keep_all = TRUE) %>% #get rid of duplicate mice since metadata encompasses samples from all timepoints
  select(row_unique_mouse, clearance_status_d7, experiment)

#Select dn1 vs d7 distances from overall distance matrix, join to mouse_clear_d7_status to add clearance_status_d7 column
dn1vd7_dist <- all_dist %>% 
  #Relabel rows and columns by separating sample id into mouse and exp labels
  separate(rows, into=c("row_mouse", "comparison_day", "row_exp"), sep = c(3, -2)) %>% #separate mouse ID from D and experiment number
  separate(columns, into=c("col_mouse", "baseline_day", "col_exp"), sep = c(3, -2)) %>% #separate mouse ID from D and experiment number
  filter(comparison_day %in% c("Dn1", "D7")) %>%  #select day -1 and day 7 timepoints
  filter(baseline_day %in% c("Dn1", "D7")) %>%  #select day -1 and day 7 timepoints
  unite("row_unique_mouse", sep = "_", row_mouse, row_exp, remove = TRUE) %>% #make columns to denote unique individual mouse from columns
  unite("col_unique_mouse", sep = "_", col_mouse, col_exp, remove = TRUE) %>% #make to denote unique individual mouse from rows
  filter(row_unique_mouse == col_unique_mouse) %>% #Select mice with samples from both D-1 and D7 timepoints. #30 mice
  select(-col_unique_mouse) %>% #get rid of since row_ and col_unique mouse labels match
  #Make a columns to designate vendors/source
  mutate(vendor = case_when(str_detect(row_unique_mouse, "S") ~ "Schloss",
                            str_detect(row_unique_mouse, "Y") ~ "Young",
                            str_detect(row_unique_mouse, "J") ~ "Jackson",
                            str_detect(row_unique_mouse, "C") ~ "Charles River",
                            str_detect(row_unique_mouse, "T") ~ "Taconic",
                            str_detect(row_unique_mouse, "E") ~ "Envigo")
  ) %>% 
  mutate(vendor=factor(vendor, levels=c("Schloss", "Young", "Jackson", "Charles River", "Taconic", "Envigo"))) %>% #Transform vendor into factor, same as what was done for metadata
  left_join(mouse_clear_d7_status, by= "row_unique_mouse")


#Plot of distances between day -1 and day 7 samples across sources/vendors
#Within experiment and same vendor comparisons
dn1_vs_d7 <- dn1vd7_dist %>% 
  group_by(vendor) %>% 
  mutate(median = median(distances)) %>% 
  ungroup() %>% 
  ggplot(aes(x = vendor, y = distances, color = vendor, shape = experiment)) +
  geom_errorbar(aes(ymax=median, ymin=median), color = "gray50", size = 1, show.legend = FALSE)+ #Add lines to indicate median 
  geom_jitter(size=2) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2))+
  scale_colour_manual(name=NULL,
                      values=color_scheme,
                      breaks=color_vendors,
                      labels=color_vendors)+
  scale_shape_manual(name=NULL,
                     values=shape_scheme,
                     breaks=shape_experiment,
                     labels=shape_experiment) +  
  labs(title = NULL, x = NULL, y = "Theta YC Distance relative to baseline")+
  theme(plot.title = element_text(hjust = 0.5)) +#Center plot title
  theme_classic()+
  theme(legend.position = "none")
save_plot("results/figures/dn1_vs_d7_thetayc.png", dn1_vs_d7)


