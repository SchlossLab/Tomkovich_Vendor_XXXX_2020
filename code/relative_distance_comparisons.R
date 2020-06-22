source("code/functions.R")

read_dist_df <- function(dist_file_name){
  linear_data <- scan(dist_file_name, what="character", quiet=TRUE)[-1]
  
  samples <- str_subset(linear_data, "D") #Pull out sample ids with E (denotes experiment number
  n_samples <- length(samples)
  distance_strings <- str_subset(linear_data, "\\.")
  
  distance_matrix <- matrix(0, nrow=n_samples, ncol=n_samples)
  colnames(distance_matrix) <- samples
  as.tibble(cbind(rows=samples, distance_matrix)) %>%
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
  geom_boxplot(outlier.shape = NA, size = 1.2)+
  geom_jitter(aes(shape = experiment), size=2, alpha=0.6, show.legend = FALSE) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2))+
  ylim(0, 1.25)+
  scale_colour_manual(name=NULL,
                      values=color_scheme,
                      breaks=color_vendors,
                      labels=color_vendors)+
  labs(title = "Baseline: within experiment", x = NULL, y = "Theta YC Distance")+
  theme(plot.title = element_text(hjust = 0.5))+ #Center plot title
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
  labs(title = "Baseline: between experiments", x = NULL, y = "Theta YC Distance")+
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
