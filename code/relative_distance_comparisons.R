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

#Within source comparisons
intra_vendors <- dn1_dist %>% 
  filter(vendor == col_vendor) %>% 
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
  labs(title = "Baseline: within group", x = NULL, y = "Theta YC Distance")+
  theme(plot.title = element_text(hjust = 0.5))+ #Center plot title
  theme_classic()
save_plot("results/figures/within_vendors_dn1.png", intra_vendors)

#Between vendor comparisons (all vendors)
inter_vendors <- dn1_dist %>% 
  filter(vendor != col_vendor) #663

#Function to select just intervendor comparison, get rid of vendor and col_vendor columns and make a new column to denote vendor
inter_vendor_select <- function(vendor_name){
  inter_vendors %>% 
  filter(vendor == vendor_name | col_vendor == vendor_name) %>% 
  select(-vendor, -col_vendor) %>% 
  mutate(vendor = vendor_name)
}
inter_schloss <- inter_vendor_select("Schloss") #204 
inter_young <- inter_vendor_select("Young") #256
inter_jackson <- inter_vendor_select("Jackson") #175
inter_charles <- inter_vendor_select("Charles River")  #204
inter_taconic <- inter_vendor_select("Taconic") #256
inter_envigo <- inter_vendor_select("Envigo") #231
#Combine intergroup comparisons for each vendor into one data frame & plot
inter_all <- rbind(inter_schloss, inter_young, inter_jackson, inter_charles,
                   inter_taconic, inter_envigo) %>% 
  mutate(vendor=factor(vendor, levels=c("Schloss", "Young", "Jackson", "Charles River", "Taconic", "Envigo"))) %>% #Transform vendor into factor, same as what was done for metadata
#  group_by(vendor) %>% 
  #  summarize(median_dist = median(distances), n = n()) %>% 
  ggplot(aes(x = vendor, y = distances, color = vendor)) +
  geom_boxplot(outlier.shape = NA, size = 1.2)+
  scale_x_discrete(guide = guide_axis(n.dodge = 2))+
  ylim(0, 1.25)+
  scale_colour_manual(name=NULL,
                      values=color_scheme,
                      breaks=color_vendors,
                      labels=color_vendors)+
  labs(title = "Baseline: between groups", x = NULL, y = "Theta YC Distance")+
  theme(plot.title = element_text(hjust = 0.5))+ #Center plot title
  theme_classic()
save_plot("results/figures/between_vendors_dn1.png", inter_all)

#Day 0 Inter and Intra-group variation----
d0_dist <- read_dist_df("data/mothur/d0/vendors.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.0.03.subsample.thetayc.0.03.lt.ave.dist") %>% 
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
within_exp_d0 <- d0_dist %>% 
  filter(vendor == col_vendor) %>% 
  filter(experiment == col_exp) %>% 
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
  labs(title = "Clindamycin: within experiment", x = NULL, y = "Theta YC Distance")+
  theme(plot.title = element_text(hjust = 0.5))+ #Center plot title
  theme_classic()
save_plot("results/figures/within_exp_d0.png", within_exp_d0)

#Between experiment and same vendor comparisons
between_exp_d0 <- d0_dist %>% 
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
  labs(title = "Clindamycin: between experiments", x = NULL, y = "Theta YC Distance")+
  theme(plot.title = element_text(hjust = 0.5))+ #Center plot title
  theme_classic()
save_plot("results/figures/between_exp_d0.png", between_exp_d0)

#Within source comparisons
intra_vendors_d0 <- d0_dist %>% 
  filter(vendor == col_vendor) %>% 
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
  labs(title = "Clindamycin: within group", x = NULL, y = "Theta YC Distance")+
  theme(plot.title = element_text(hjust = 0.5))+ #Center plot title
  theme_classic()
save_plot("results/figures/within_vendors_d0.png", intra_vendors_d0)

#Between vendor comparisons (all vendors)
inter_vendors_d0 <- d0_dist %>% 
  filter(vendor != col_vendor) #629

#Function to select just intervendor comparison, get rid of vendor and col_vendor columns and make a new column to denote vendor
inter_vendor_select <- function(vendor_name){
  inter_vendors_d0 %>% 
    filter(vendor == vendor_name | col_vendor == vendor_name) %>% 
    select(-vendor, -col_vendor) %>% 
    mutate(vendor = vendor_name)
}
inter_schloss_d0 <- inter_vendor_select("Schloss")  
inter_young_d0 <- inter_vendor_select("Young") 
inter_jackson_d0 <- inter_vendor_select("Jackson") 
inter_charles_d0 <- inter_vendor_select("Charles River") 
inter_taconic_d0 <- inter_vendor_select("Taconic") 
inter_envigo_d0 <- inter_vendor_select("Envigo") 
#Combine intergroup comparisons for each vendor into one data frame & plot
inter_all_d0 <- rbind(inter_schloss_d0, inter_young_d0, inter_jackson_d0, inter_charles_d0,
                   inter_taconic_d0, inter_envigo_d0) %>% 
  mutate(vendor=factor(vendor, levels=c("Schloss", "Young", "Jackson", "Charles River", "Taconic", "Envigo"))) %>% #Transform vendor into factor, same as what was done for metadata
  #  group_by(vendor) %>% 
  #  summarize(median_dist = median(distances), n = n()) %>% 
  ggplot(aes(x = vendor, y = distances, color = vendor)) +
  geom_boxplot(outlier.shape = NA, size = 1.2)+
  scale_x_discrete(guide = guide_axis(n.dodge = 2))+
  ylim(0, 1.25)+
  scale_colour_manual(name=NULL,
                      values=color_scheme,
                      breaks=color_vendors,
                      labels=color_vendors)+
  labs(title = "Clindamycin: between groups", x = NULL, y = "Theta YC Distance")+
  theme(plot.title = element_text(hjust = 0.5))+ #Center plot title
  theme_classic()  
save_plot("results/figures/between_vendors_d0.png", inter_all_d0)

#Day 1 Inter and Intra-group variation----
d1_dist <- read_dist_df("data/mothur/d1/vendors.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.0.03.subsample.thetayc.0.03.lt.ave.dist") %>% 
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
within_exp_d1 <- d1_dist %>% 
  filter(vendor == col_vendor) %>% 
  filter(experiment == col_exp) %>% 
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
  labs(title = "Post-infection: within experiment", x = NULL, y = "Theta YC Distance")+
  theme(plot.title = element_text(hjust = 0.5))+ #Center plot title
  theme_classic()
save_plot("results/figures/within_exp_d1.png", within_exp_d1)

#Between experiment and same vendor comparisons
between_exp_d1 <- d1_dist %>% 
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
  labs(title = "Post-infection: between experiments", x = NULL, y = "Theta YC Distance")+
  theme(plot.title = element_text(hjust = 0.5))+ #Center plot title
  theme_classic()
save_plot("results/figures/between_exp_d1.png", between_exp_d1)

#Within source comparisons
intra_vendors_d1 <- d1_dist %>% 
  filter(vendor == col_vendor) %>% 
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
  labs(title = "Post-infection: within group", x = NULL, y = "Theta YC Distance")+
  theme(plot.title = element_text(hjust = 0.5))+ #Center plot title
  theme_classic()
save_plot("results/figures/within_vendors_d1.png", intra_vendors_d1)

#Between vendor comparisons (all vendors)
inter_vendors_d1 <- d1_dist %>% 
  filter(vendor != col_vendor) #629

#Function to select just intervendor comparison, get rid of vendor and col_vendor columns and make a new column to denote vendor
inter_vendor_select <- function(vendor_name){
  inter_vendors_d1 %>% 
    filter(vendor == vendor_name | col_vendor == vendor_name) %>% 
    select(-vendor, -col_vendor) %>% 
    mutate(vendor = vendor_name)
}
inter_schloss_d1 <- inter_vendor_select("Schloss") 
inter_young_d1 <- inter_vendor_select("Young") 
inter_jackson_d1 <- inter_vendor_select("Jackson") 
inter_charles_d1 <- inter_vendor_select("Charles River")  
inter_taconic_d1 <- inter_vendor_select("Taconic") 
inter_envigo_d1 <- inter_vendor_select("Envigo") 
#Combine intergroup comparisons for each vendor into one data frame & plot
inter_all_d1 <- rbind(inter_schloss_d1, inter_young_d1, inter_jackson_d1, inter_charles_d1,
                      inter_taconic_d1, inter_envigo_d1) %>% 
  mutate(vendor=factor(vendor, levels=c("Schloss", "Young", "Jackson", "Charles River", "Taconic", "Envigo"))) %>% #Transform vendor into factor, same as what was done for metadata
  #  group_by(vendor) %>% 
  #  summarize(median_dist = median(distances), n = n()) %>% 
  ggplot(aes(x = vendor, y = distances, color = vendor)) +
  geom_boxplot(outlier.shape = NA, size = 1.2)+
  scale_x_discrete(guide = guide_axis(n.dodge = 2))+
  ylim(0, 1.25)+
  scale_colour_manual(name=NULL,
                      values=color_scheme,
                      breaks=color_vendors,
                      labels=color_vendors)+
  labs(title = "Clindamycin: between groups", x = NULL, y = "Theta YC Distance")+
  theme(plot.title = element_text(hjust = 0.5))+ #Center plot title
  theme_classic()
save_plot("results/figures/between_vendors_d1.png", inter_all_d1)

#Day -1 versus D7 thetayc distances in individual mice, compared across sources----
#read in theta YC distances for all samples, all timepoints
all_dist <- read_dist_df("data/process/vendors.subsample.thetayc.ave.dist")

#Modify metadata to match row_unique_mouse, select clearance_status_d7 column and join to dn1vd7_dist
mouse_clear_d7_status <- metadata %>% 
  separate(id, into=c("row_mouse", "comparison_day", "row_exp"), sep = c(3, -2)) %>% #separate mouse ID from D and experiment number
  unite("row_unique_mouse", sep = "_", row_mouse, row_exp, remove = TRUE) %>% #make columns to denote unique individual mouse from columns
  distinct(row_unique_mouse, .keep_all = TRUE) %>% #get rid of duplicate mice since metadata encompasses samples from all timepoints
  select(row_unique_mouse, clearance_status_d7)

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
  ggplot(aes(x = vendor, y = distances, color = vendor)) +
  geom_errorbar(aes(ymax=median, ymin=median), color = "gray50", size = 1, show.legend = FALSE)+ #Add lines to indicate median 
  geom_jitter(size=2) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2))+
  scale_colour_manual(name=NULL,
                      values=color_scheme,
                      breaks=color_vendors,
                      labels=color_vendors)+
#Scale for adding aes(shape = clearance_status_d7) to geom_jitter...  
#  scale_shape_manual(name="Cleared by Day 7",
#                     values=c(4, 19, 21),
#                     breaks=c("colonized", "not_detectable", "no_data"),
#                     labels=c("no", "yes", "no data"), 
#                     drop=FALSE, na.translate = TRUE, na.value = 1)+
  labs(title = NULL, x = NULL, y = "Theta YC Distance relative to baseline")+
  theme(plot.title = element_text(hjust = 0.5)) +#Center plot title
  theme_classic()+
  theme(axis.text.x = element_blank()) #Remove X-axis labels because they are redundant with our key
save_plot("results/figures/dn1_vs_d7_thetayc.png", dn1_vs_d7)

#Theta YC distances over time relative to baseline (dn1)

#Function to examine distances at a specific timepoint relative to baseline
#Argument: timepoint to examine distances relative to baseline
day_spec_dist <- function(timepoint){
  all_dist %>% 
    #Relabel rows and columns by separating sample id into mouse and exp labels
    separate(rows, into=c("row_mouse", "comparison_day", "row_exp"), sep = c(3, -2)) %>% #separate mouse ID from D and experiment number
    separate(columns, into=c("col_mouse", "baseline_day", "col_exp"), sep = c(3, -2)) %>% #separate mouse ID from D and experiment number
    filter(comparison_day %in% c("Dn1", timepoint)) %>%  #select day -1 and day 7 timepoints
    filter(baseline_day %in% c("Dn1", timepoint)) %>%  #select day -1 and day 7 timepoints
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
    mutate(vendor=factor(vendor, levels=c("Schloss", "Young", "Jackson", "Charles River", "Taconic", "Envigo"))) #Transform vendor into factor, same as what was done for metadata
}

d0_dist <-day_spec_dist("D0")  
d1_dist <-day_spec_dist("D1")  
d2_dist <-day_spec_dist("D2")  
d3_dist <-day_spec_dist("D3") 
d4_dist <-day_spec_dist("D4") 
d5_dist <-day_spec_dist("D5") 
d6_dist <-day_spec_dist("D6") 
d7_dist <-day_spec_dist("D7") 
d8_dist <-day_spec_dist("D8") 
d9_dist <-day_spec_dist("D9")

rel_dist_over_time <- rbind(d0_dist, d1_dist, d2_dist, d3_dist, d4_dist, d5_dist, d6_dist, d7_dist, d8_dist, d9_dist) %>% 
  mutate(comparison_day = recode(comparison_day, D0 = 0, D1 = 1, D2 = 2, D3 = 3,
                                 D4 = 4, D5 = 5, D6 = 6, D7 = 7, D8 = 8, D9 = 9))  #Match experiment column variable names from metadata
  #fix comparison day column to get rid of D before timepoint

#Number of mice with paired samples per timepoint:
n_per_timepoint <- rel_dist_over_time %>% 
  group_by(comparison_day) %>% count()

#Plot of how thetayc distances relative to baseline change over time across each source
#Note need to get clearance status day 7 for each mouse from metadata first
rel_dist_median <- rel_dist_over_time %>% 
  group_by(vendor, comparison_day) %>% 
  summarize(median=(median(distances))) %>% 
  ungroup
plot_rel_dist_over_time <-  ggplot(NULL)+
    geom_point(rel_dist_over_time, mapping = aes(x=comparison_day, y=distances, color=vendor, shape=clearance_status_d7), size  = 1.5, position = position_dodge(width = 0.6))+
    geom_line(rel_dist_median, mapping = aes(x=comparison_day, y=median, color=vendor), size = 1, show.legend = FALSE)+
    scale_colour_manual(name=NULL,
                        values=color_scheme,
                        breaks=color_vendors,
                        labels=color_vendors)+
    scale_shape_manual(name="Cleared by Day 7",
                       values=c(4, 19, 21),
                       breaks=c("colonized", "not_detectable", "no_data"),
                       labels=c("no", "yes", "no data"), 
                       drop=FALSE, na.translate = TRUE, na.value = 1)+
    labs(title=NULL,
         x="Day",
         y="Theta YC Distance relative to baseline") +
    scale_x_continuous(breaks = c(-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9),
                       limits = c(-1.5, 9.5)) +
    theme_classic()+
    theme(plot.title=element_markdown(hjust = 0.5),
          text = element_text(size = 16)) # Change font size for entire plot

#Plot without clearance status day 7 shape aesthetic
rel_dist_median <- rel_dist_over_time %>% 
  group_by(vendor, comparison_day) %>% 
  summarize(median=(median(distances))) %>% 
  ungroup
plot_rel_dist_over_time <-  ggplot(NULL)+
  geom_point(rel_dist_over_time, mapping = aes(x=comparison_day, y=distances, color=vendor), size  = 1.5, position = position_dodge(width = 0.6))+
  geom_line(rel_dist_median, mapping = aes(x=comparison_day, y=median, color=vendor), size = 1, show.legend = FALSE)+
  scale_colour_manual(name=NULL,
                      values=color_scheme,
                      breaks=color_vendors,
                      labels=color_vendors)+
  labs(title=NULL,
       x="Day",
       y="Theta YC Distance relative to baseline") +
  scale_x_continuous(breaks = c(-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9),
                     limits = c(-1.5, 9.5)) +
  theme_classic()+
  theme(plot.title=element_markdown(hjust = 0.5),
        text = element_text(size = 16)) # Change font size for entire plot


