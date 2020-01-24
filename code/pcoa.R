source("code/functions.R")

pcoa_data <- read_tsv("data/process/vendors.subsample.thetayc.ave.pcoa.axes") %>%
  select(group, axis1, axis2) %>% #Limit to 2 PCoA axes
  rename(id = group) %>% #group is the same as id in the metadata data frame
  right_join(metadata, by= "id") #merge metadata and PCoA data frames
  
#Function to plot pcoa data for all vendors----
plot_pcoa <- function(df){
  ggplot(df, aes(x=axis1, y=axis2, color = vendor, alpha = day)) +
	geom_point(size=2) +
  scale_colour_manual(name=NULL,
                      values=color_scheme,
                      breaks=color_vendors,
                      labels=color_vendors)+
  scale_alpha_continuous(range = c(.3, 1),
                         breaks= c(2, 4, 6, 8, 10),
                         labels=c(2, 4, 6, 8, 10))+
	coord_fixed() + 
	labs(x="PCoA 1",
	     y="PCoA 2",
	     color= "Vendor",
	     alpha= "Day") +
	theme_classic()
}

#PCoA plot that combines the 2 experiments----  
pcoa_plot_combined <- plot_pcoa(pcoa_data)

#PCoA plot for the 1st experiment----
pcoa_exp1 <- plot_pcoa(pcoa_data %>% filter(experiment == 1))

#PCoA plot for the 2nd experiment----
pcoa_exp2 <- plot_pcoa(pcoa_data %>% filter(experiment == 2))

plot_grid(pcoa_plot_combined, pcoa_exp1, pcoa_exp2, labels = c("Combined Experiments", "Experiment 1", "Experiment 2"), ncol =1, label_x = .086, label_y = 1)+
  ggsave("exploratory/notebook/pcoa_by_exp.pdf", width = 5, height = 11)

#Function to plot pcoa data for all vendors----
plot_pcoa_timepoint <- function(df, desired_day){
  pcoa_timeppoint <- ggplot(df, aes(x=axis1, y=axis2, color = vendor)) +
    geom_point(size=2, alpha = 0.4) +
    scale_colour_manual(name=NULL,
                        values=color_scheme,
                        breaks=color_vendors,
                        labels=color_vendors)+
    coord_fixed() + 
    xlim(-0.4, 0.6)+
    ylim(-0.6, 0.4)+
    labs(x="PCoA 1",
         y="PCoA 2",
         color= "Vendor",
         alpha= "Day") +
    theme_classic()
  save_plot(filename = paste0("results/figures/pcoa_day", desired_day,".png"), pcoa_timeppoint)
}
#PCoA plot for D0, 1, 2, 4, 5, 8, 9 timepoints----
pcoa_day0 <- plot_pcoa_timepoint(pcoa_data %>% filter(day == 0), 0)
pcoa_day1 <- plot_pcoa_timepoint(pcoa_data %>% filter(day == 1), 1)
pcoa_day2 <- plot_pcoa_timepoint(pcoa_data %>% filter(day == 2), 2)
pcoa_day4 <- plot_pcoa_timepoint(pcoa_data %>% filter(day == 4), 4)
pcoa_day5 <- plot_pcoa_timepoint(pcoa_data %>% filter(day == 5), 5)  
pcoa_day8 <- plot_pcoa_timepoint(pcoa_data %>% filter(day == 8), 8)
pcoa_day9 <- plot_pcoa_timepoint(pcoa_data %>% filter(day == 9), 9)

#Plot pcoa data by experiment on Day -1 (start of the experiment before clindamycin treatment) ----
plot_initial_communities <- pcoa_data %>% filter(day == -1) %>% 
  ggplot(aes(x=axis1, y=axis2, color = vendor, shape = experiment)) +
  scale_colour_manual(name=NULL,
                      values=color_scheme,
                      breaks=color_vendors,
                      labels=color_vendors)+
    geom_point(size=2) +
    coord_fixed() + 
    labs(x="PCoA 1",
         y="PCoA 2",
         color= "Vendor",
         shape= "Experiment") +
    theme_classic()+
  ggsave("exploratory/notebook/pcoa_initial_all.pdf")

#Plot each vendor separately for the 2 experiments (all timepoints)----

#Function to plot pcoa data for each vendor separately----
plot_pcoa_vendor <- function(indiv_vendor){
  pcoa_data %>% filter(vendor == indiv_vendor) %>% 
  ggplot(aes(x=axis1, y=axis2, color = experiment, alpha = day)) +
    geom_point(size=2) +
    scale_alpha_continuous(range = c(.3, 1),
                           breaks= c(2, 4, 6, 8, 10),
                           labels=c(2, 4, 6, 8, 10))+
    coord_fixed() + 
    xlim(-0.4, 0.6)+
    ylim(-0.6, 0.4)+
    labs(title = indiv_vendor,
         x="PCoA 1",
         y="PCoA 2",
         color= "Experiment",
         alpha= "Day") +
    theme_classic()
}

Schloss <- plot_pcoa_vendor("Schloss")
Young <- plot_pcoa_vendor("Young")
Charles_River <- plot_pcoa_vendor("Charles River") 
Envigo <- plot_pcoa_vendor("Envigo")
Jackson <- plot_pcoa_vendor("Jackson")
Taconic <- plot_pcoa_vendor("Taconic")

plot_grid(Schloss, Young, Charles_River, Envigo, Jackson, Taconic,
          ncol =2, label_x = .086, label_y = 1)+
  ggsave("exploratory/notebook/pcoa_by_vendor.pdf", width = 8.5, height = 11)

#Function to plot pcoa data for each vendor separately for initial timepoint only----
plot_pcoa_vendor_initial <- function(indiv_vendor){
  pcoa_data %>% filter(vendor == indiv_vendor) %>% 
    filter(day == -1) %>% 
    ggplot(aes(x=axis1, y=axis2, color = experiment)) +
    geom_point(size=2) +
    coord_fixed() + 
    xlim(-0.4, 0.6)+
    ylim(-0.6, 0.4)+
    labs(title = indiv_vendor,
         x="PCoA 1",
         y="PCoA 2",
         color= "Experiment",
         alpha= "Day") +
    theme_classic()
}

Schloss_initial <- plot_pcoa_vendor_initial("Schloss")
Young_initial <- plot_pcoa_vendor_initial("Young")
Charles_River_initial <- plot_pcoa_vendor_initial("Charles River") 
Envigo_initial <- plot_pcoa_vendor_initial("Envigo")
Jackson_initial <- plot_pcoa_vendor_initial("Jackson")
Taconic_initial <- plot_pcoa_vendor_initial("Taconic")

plot_grid(Schloss_initial, Young_initial, Charles_River_initial, Envigo_initial, Jackson_initial, Taconic_initial,
          ncol = 2, label_x = .086, label_y = 1)+
  ggsave("exploratory/notebook/pcoa_by_vendor_initial.pdf", width = 8.5, height = 11)

#Check duplicate samples on PCoA----
pcoa_data_duplicates <- read_tsv("data/process/vendors.subsample.thetayc.ave.pcoa.axes") %>%
  select(group, axis1, axis2) %>% #Limit to 2 PCoA axes
  rename(id = group) %>% #group is the same as id in the metadata data frame
  left_join(metadata, by= "id") #merge metadata and PCoA data frames

plot_pcoa_duplicates <- function(limited_dataframe){
    limited_dataframe %>% 
    ggplot(aes(x=axis1, y=axis2, label = id)) +
    geom_point(size=3, , alpha = 0.4) + 
    geom_text(color = "black", size =2,
              hjust = 0, nudge_x = 0.05)+
    coord_fixed() + 
    xlim(-0.4, 0.6)+
    ylim(-0.6, 0.4)+
    labs(title = NULL,
         x="PCoA 1",
         y="PCoA 2") +
    theme_classic()
}

#Check list of duplicates identified in functions.R. Presented these in lab meeting
C21D0E2 <- plot_pcoa_duplicates(pcoa_data_duplicates %>% filter(id == 'C21D0E2'| id == 'C21D0E2no2'))
E212Dn1E2 <- plot_pcoa_duplicates(pcoa_data_duplicates %>% filter(id == 'E212Dn1E2'| id == 'E21Dn1E2'))
E21D1E2 <- plot_pcoa_duplicates(pcoa_data_duplicates %>% filter(id == 'E21D1E2No1'| id == 'E21D1E2No2'))
E21Dn1E1no2 <- plot_pcoa_duplicates(pcoa_data_duplicates %>% filter(id == 'E21Dn1E1no2'| id == 'E21Dn1E1'))
S22D1E2 <- plot_pcoa_duplicates(pcoa_data_duplicates %>% filter(id == 'S22D1E2No1'| id == 'S22D1E2No2'))
T12D8E1 <- plot_pcoa_duplicates(pcoa_data_duplicates %>% filter(id == 'T12D8E1No1'| id == 'T12D8E1No2'))
Y13D8E2 <- plot_pcoa_duplicates(pcoa_data_duplicates %>% filter(id == 'Y13D8E2'| id == 'Y13D9E1'))

plot_grid(C21D0E2, E212Dn1E2, E21D1E2, E21Dn1E1no2, S22D1E2, T12D8E1, Y13D8E2, labels = c("C21D0E2", "E212Dn1E2", "E21D1E2", "E21Dn1E1no2", "S22D1E2", "T12D8E1", "Y13D8E2"),
          ncol = 2, label_x = .4, label_y = 1)+
  ggsave("exploratory/notebook/pcoa_duplicates.pdf", width = 8.5, height = 11)

#Plot duplicates in context of the rest of the mice from same group & timepoint----
plot_pcoa_duplicates_context <- function(limited_dataframe){
  limited_dataframe %>% 
    ggplot(aes(x=axis1, y=axis2, label = id, color = experiment)) +
    geom_point(size=3, , alpha = 0.4) + 
    geom_text(color = "black", size =2,
              hjust = 0, nudge_x = 0.05)+
    coord_fixed() + 
    xlim(-0.4, 0.6)+
    ylim(-0.6, 0.4)+
    labs(title = NULL,
         x="PCoA 1",
         y="PCoA 2") +
    theme_classic()
}

## Plots of duplicate samples in the context of the rest of the samples from that group & timepoint----
# C21D0E2no2 (duplicate of C21D0E2?, which is present in shared_sample_names (both part of run_1, plate_3) or C21D0E1 which isn't present)
C21D0 <- plot_pcoa_duplicates_context(pcoa_data_duplicates %>% filter(id == 'C21D0E2'| id == 'C21D0E2no2' | vendor == "Charles River" & day == "0")) 
# 6 samples total

# E212Dn1E2 (duplicate of E21Dn1E2?, which is present in shared_sample_names or E12Dn1E2, which is absent in shared_sample_names? E22Dn1E2, E12Dn1E1, E11Dn1E2 is also present in shared_sample names).
E212Dn1 <- plot_pcoa_duplicates_context(pcoa_data_duplicates %>% filter(id == 'E212Dn1E2'| id == 'E21Dn1E2' | vendor == "Envigo" & day == "-1"))

# E21D1E2No1 and E21D1E2No2. Sample sequenced twice? Both were part of the same plate on run_1, plate_2.
E21D1 <- plot_pcoa_duplicates_context(pcoa_data_duplicates %>% filter(id == 'E21D1E2No1'| id == 'E21D1E2No2' | vendor == "Envigo" & day == "1"))

# E21Dn1E1no2 (duplicate of E21Dn1E1?, which is present in shared_sample_names). Both a part of run_1, plate_3 and E21Dn1E2 is also present on run_1, plate_3
E21Dn1 <- plot_pcoa_duplicates_context(pcoa_data_duplicates %>% filter(id == 'E21Dn1E1no2'| id == 'E21Dn1E1' | vendor == "Envigo" & day == "-1"))

# S22D1E2No1 and S22D1E2No2. Sample sequenced twice? Which one to pick...Both a part of run_2, plate_1 and S22D1E2 was also sequenced on run_2, plate_3
S22D1 <- plot_pcoa_duplicates_context(pcoa_data_duplicates %>% filter(id == 'S22D1E2No1'| id == 'S22D1E2No2' | vendor == "Schloss" & day == "1"))

# T12D8E1No1 and T12D8E1No2. Sample sequenced twice? Which one to pick...Both  a part of run_2, plate_3, T12D8E2 is in shared_sample_names but not run_plate_info? 
T12D8 <- plot_pcoa_duplicates_context(pcoa_data_duplicates %>% filter(id == 'T12D8E1No1'| id == 'T12D8E1No2' | vendor == "Taconic" & day == "8"))

# Y13D8E2. Okay to lose this sample? There was only a Y13 mouse in the 1st experiment and Y13D8E1 is already listed in shared_sample_names. 
#Looking back at the plate layout maps Y13D9E1 (Y13_D9_E1) was listed twice, so maybe this sample was really a duplicate of that one.
#Y13D9E1 shows up in run_plate_info duplicates, while Y13D8E2 is not present in run_plate_info
Y13D8 <- plot_pcoa_duplicates_context(pcoa_data_duplicates %>% filter(id == 'Y13D8E2'| id == 'Y13D9E1' | vendor == "Young" & day == "9"))

plot_grid(C21D0, E212Dn1, E21D1, E21Dn1, S22D1, T12D8, Y13D8, labels = c("CR Day 0", "Envigo Day -1", "Envigo Day 1", "Envigo Day -1", "Schloss Day 1", "Taconic D8", "Young Day 9"),
          ncol = 2, label_x = .4, label_y = 1)+
  ggsave("exploratory/notebook/pcoa_duplicates_in_context.pdf", width = 8.5, height = 11)

