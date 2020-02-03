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
    ylim(-0.6, 0.41)+
    labs(x="PCoA 1",
         y="PCoA 2",
         color= "Vendor",
         alpha= "Day") +
    theme_classic()
  save_plot(filename = paste0("results/figures/pcoa_day", desired_day,".png"), pcoa_timeppoint)
}
#PCoA plot for D0, 1, 2, 3, 4, 5, 6, 7, 8, 9 timepoints----
pcoa_dayn1 <- plot_pcoa_timepoint(pcoa_data %>% filter(day == -1), -1)
pcoa_day0 <- plot_pcoa_timepoint(pcoa_data %>% filter(day == 0), 0)
pcoa_day1 <- plot_pcoa_timepoint(pcoa_data %>% filter(day == 1), 1)
pcoa_day2 <- plot_pcoa_timepoint(pcoa_data %>% filter(day == 2), 2)
pcoa_day3 <- plot_pcoa_timepoint(pcoa_data %>% filter(day == 3), 3)
pcoa_day4 <- plot_pcoa_timepoint(pcoa_data %>% filter(day == 4), 4)
pcoa_day5 <- plot_pcoa_timepoint(pcoa_data %>% filter(day == 5), 5)
pcoa_day6 <- plot_pcoa_timepoint(pcoa_data %>% filter(day == 6), 6)
pcoa_day7 <- plot_pcoa_timepoint(pcoa_data %>% filter(day == 7), 7)
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
    ylim(-0.6, 0.41)+
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
    ylim(-0.6, 0.41)+
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

#Check if points separate according to MiSeq run number.----
#Added MiSeq run #, plate# and plate_location columns to metadata in functions.R script. (Created 23 duplicate sample rows because these were sequenced twice on separate runs)
plot_pcoa_miseq_run <- pcoa_data %>% 
    ggplot(aes(x=axis1, y=axis2, color = run)) +
    geom_point(size=2, alpha = 0.4) +
#    scale_alpha_continuous(range = c(.3, 1),
#                           breaks= c(2, 4, 6, 8, 10),
#                           labels=c(2, 4, 6, 8, 10))+
    coord_fixed() + 
    xlim(-0.4, 0.6)+
    ylim(-0.6, 0.41)+
    labs(title = NULL,
         x="PCoA 1",
         y="PCoA 2",
         color= "MiSeq Run",
         alpha= "Day") +
    theme_classic()+
  ggsave("exploratory/notebook/pcoa_by_MiSeq_run.pdf")
#All the NA samples for run are from D8 of Experiment 2

#Check within run duplicate samples on PCoA----
pcoa_data_duplicates <- read_tsv("data/process/vendors.subsample.thetayc.ave.pcoa.axes") %>%
  select(group, axis1, axis2) %>% #Limit to 2 PCoA axes
  rename(id = group) %>% #group is the same as id in the metadata data frame
  left_join(metadata, by= "id") #merge metadata and PCoA data frames

plot_pcoa_duplicates <- function(limited_dataframe){
    limited_dataframe %>% 
    ggplot(aes(x=axis1, y=axis2, label = id)) +
    geom_point(size=3, alpha = 0.2) + 
    geom_text(color = "black", size =2,
              hjust = 0, nudge_x = 0.05)+
    coord_fixed() + 
    xlim(-0.4, 0.6)+
    ylim(-0.6, 0.41)+
    labs(title = NULL,
         x="PCoA 1",
         y="PCoA 2") +
    theme_classic()
}

#Check list of within run duplicates identified in functions.R. Presented these in lab meeting
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

#Plot within run duplicates in context of the rest of the mice from same group & timepoint----
plot_pcoa_duplicates_context <- function(limited_dataframe){
  limited_dataframe %>% 
    ggplot(aes(x=axis1, y=axis2, label = id, color = experiment)) +
    geom_point(size=3, alpha = 0.4) + 
    geom_text(color = "black", size =2,
              hjust = 0, nudge_x = 0.05)+
    coord_fixed() + 
    xlim(-0.4, 0.6)+
    ylim(-0.6, 0.41)+
    labs(title = NULL,
         x="PCoA 1",
         y="PCoA 2") +
    theme_classic()
}

## Plots of within run duplicate samples in the context of the rest of the samples from that group & timepoint----
# C21D0E2no2 (duplicate of C21D0E2?, which is present in shared_sample_names (both part of run_1, plate_3) or C21D0E1 which isn't present)
C21D0 <- plot_pcoa_duplicates_context(pcoa_data_duplicates %>% filter(id == 'C21D0E2'| id == 'C21D0E2no2' | vendor == "Charles River" & day == "0")) # 6 samples total
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

#Examine the sample ids that were missing from plate_map but run through mothur on a PCoA----
missing_platemap <- c("C11D8E2", "C21D8E2", "C22D8E2", "E12D8E2", "E21D8E2", "E22D8E2", "J11D8E2", "J21D8E2", "J22D8E2", "S12D8E2", "S21D8E2", "S22D8E2", "T11D8E2", "T12D8E2", "T21D8E2", "Y12D8E2", "Y13D8E2", "Y21D8E2", "Y22D8E2") 
#All have corresponding .fastq files in the raw data folder but not a corresponding pair of .fastq files in the BaseCalls folder of either sequencing run
plot_pcoa_missing_platemap <-   function(df){
  df %>% 
  ggplot(aes(x=axis1, y=axis2, color = vendor)) +
    geom_point(size=2) +
    scale_colour_manual(name=NULL,
                        values=color_scheme,
                        breaks=color_vendors,
                        labels=color_vendors)+
    coord_fixed() + 
    xlim(-0.4, 0.6)+
    ylim(-0.6, 0.41)+
    labs(x="PCoA 1",
         y="PCoA 2",
         color= "Vendor") +
    theme_classic()
}

missing_platemap_pcoa <- plot_pcoa_missing_platemap(pcoa_data %>% filter(id %in% missing_platemap))  
plot_pcoa_d7E1 <- plot_pcoa_missing_platemap(pcoa_data %>% filter(day == 7 & experiment == 1)) 
plot_pcoa_d7E2 <- plot_pcoa_missing_platemap(pcoa_data %>% filter(day == 7 & experiment == 2)) 
plot_pcoa_d8E1 <- plot_pcoa_missing_platemap(pcoa_data %>% filter(day == 8 & experiment == 1)) 
plot_pcoa_d8E2 <- plot_pcoa_missing_platemap(pcoa_data %>% filter(day == 8 & experiment == 2))
plot_pcoa_d9E1 <- plot_pcoa_missing_platemap(pcoa_data %>% filter(day == 9 & experiment == 1))
plot_pcoa_d9E2 <- plot_pcoa_missing_platemap(pcoa_data %>% filter(day == 9 & experiment == 2))

#Examine duplicates that were sequenced on separate runs on PCoA----
duplicates_across_runs <-  c("J12D9E1", "S11D9E1", "C12D9E1", "J22D9E1", "S21D9E1", "J11D9E1", "C21D9E1", "T11D9E1", "S22D9E1", "E22D9E1", "E21D9E1", "Y13D9E1", "T21D9E1", "S12D9E1", "J21D9E1", "C11D9E1", "Y12D9E1", "T12D9E1", "Y21D9E1", "E12D9E1", "C22D9E1", "E11D9E1", "Y22D9E1")
duplicates_across_runs_pcoa <- plot_pcoa_missing_platemap(pcoa_data %>% filter(id %in% duplicates_across_runs))

plot_grid(missing_platemap_pcoa, duplicates_across_runs_pcoa, plot_pcoa_d7E1, plot_pcoa_d7E2, plot_pcoa_d8E1, plot_pcoa_d8E2, plot_pcoa_d9E1, plot_pcoa_d9E2, 
          labels = c("Samples not on MiSeq Plate Maps\n all Day 8 Experiment 2", "23 Duplicates Across MiSeq Runs\n all Day 9 Experiment 1", "Day 7 Experiment 1", "Day 7 Experiment 2", "Day 8 Experiment 1", "Day 8 Experiment 2", "Day 9 Experiment 1", "Day 9 Experiment 2"),
          ncol = 2, label_x = 0, label_y = 1)+
  ggsave("exploratory/notebook/pcoa_missing_from platemap & duplicates_across_runs.pdf", width = 8.5, height = 11)

#Directly compare Day 8 experiment 2 samples (19 total, 4 lost during subsampling) to Day 9 experiment 1 samples that were duplicates (23 total)----
# The 4 Day 8 experiment 2 samples that don't have a corresponding match were dropped during subsampling ("C12D8E2" "E11D8E2" "J12D8E2" "S11D8E2").
C11D8E2vD9E1 <- plot_pcoa_duplicates(pcoa_data_duplicates %>% filter(id == 'C11D8E2'| id == 'C11D9E1'))
C21D8E2vD9E1 <- plot_pcoa_duplicates(pcoa_data_duplicates %>% filter(id == 'C21D8E2'| id == 'C21D9E1'))
C22D8E2vD9E1 <- plot_pcoa_duplicates(pcoa_data_duplicates %>% filter(id == 'C22D8E2'| id == 'C22D9E1'))
#E11D8E2 dropped during rarefaction, but should not exist because it was found dead prior to D3
E12D8E2vD9E1 <-  plot_pcoa_duplicates(pcoa_data_duplicates %>% filter(id == 'E12D8E2'| id == 'E12D9E1'))
E21D8E2vD9E1 <-  plot_pcoa_duplicates(pcoa_data_duplicates %>% filter(id == 'E21D8E2'| id == 'E21D9E1'))
E22D8E2vD9E1 <- plot_pcoa_duplicates(pcoa_data_duplicates %>% filter(id == 'E22D8E2'| id == 'E22D9E1'))
J11D8E2vD9E1 <- plot_pcoa_duplicates(pcoa_data_duplicates %>% filter(id == 'J11D8E2'| id == 'J11D9E1'))
#J21D8E2 should not exist because it was found dead prior to Day 3
J21D8E2vD9E1 <- plot_pcoa_duplicates(pcoa_data_duplicates %>% filter(id == 'J21D8E2'| id == 'J21D9E1'))
J22D8E2vD9E1 <- plot_pcoa_duplicates(pcoa_data_duplicates %>% filter(id == 'J22D8E2' | id == 'J22D9E1'))
S12D8E2vD9E1 <- plot_pcoa_duplicates(pcoa_data_duplicates %>% filter(id == 'S12D8E2' | id == 'S12D9E1'))
S21D8E2vD9E1 <- plot_pcoa_duplicates(pcoa_data_duplicates %>% filter(id == 'S21D8E2' | id == 'S21D9E1'))
S22D8E2vD9E1 <- plot_pcoa_duplicates(pcoa_data_duplicates %>% filter(id == 'S22D8E2' | id == 'S22D9E1'))
T11D8E2vD9E1 <- plot_pcoa_duplicates(pcoa_data_duplicates %>% filter(id == 'T11D8E2' | id == 'T11D9E1'))
T12D8E2vD9E1 <- plot_pcoa_duplicates(pcoa_data_duplicates %>% filter(id == 'T12D8E2' | id == 'T12D9E1'))
T21D8E2vD9E1 <- plot_pcoa_duplicates(pcoa_data_duplicates %>% filter(id == 'T21D8E2' | id == 'T21D9E1'))
#T22 was NA in Experiment 1 but should have been present in Experiment 2
#Y11 is missing for D8 Experiment 2
Y12D8E2vD9E1 <- plot_pcoa_duplicates(pcoa_data_duplicates %>% filter(id == 'Y12D8E2' | id == 'Y12D9E1'))
Y13D8E2vD9E1 <- plot_pcoa_duplicates(pcoa_data_duplicates %>% filter(id == 'Y13D8E2' | id == 'Y13D9E1'))
#Y13 should only exist in Experiment 1
Y21D8E2vD9E1 <- plot_pcoa_duplicates(pcoa_data_duplicates %>% filter(id == 'Y21D8E2' | id == 'Y21D9E1'))
Y22D8E2vD9E1 <- plot_pcoa_duplicates(pcoa_data_duplicates %>% filter(id == 'Y22D8E2' | id == 'Y22D9E1'))

#Compare trajectories of individual mice within set of duplicate samples assuming they are part of experiment 1----
#Function to plot duplicates across runs in context of that same mouse over the rest of the timepoints from experiment 1
pcoa_duplicates_time_context <- function(limited_dataframe){
  limited_dataframe %>% 
    ggplot(aes(x=axis1, y=axis2, label = id, alpha = day, color = experiment)) +
    geom_point(size=3, show.legend = FALSE) + 
    scale_alpha_continuous(range = c(.3, 1),
                           breaks= c(2, 4, 6, 8, 10),
                           labels=c(2, 4, 6, 8, 10),
                           guide = "none")+
    geom_text(color = "black", size =2,
              hjust = 0.5, vjust = 0.5, nudge_x = 0.034)+
    coord_fixed() + 
#    xlim(-0.4, 0.6)+
#    ylim(-0.6, 0.41)+
    labs(title = NULL,
         x="PCoA 1",
         y="PCoA 2") +
    theme_classic()
}

#C11 Mouse----
C11_E1 <- pcoa_duplicates_time_context(pcoa_data_duplicates %>% filter(id == 'C11D8E2'| id == 'C11D9E1' | mouse_id == "1_Charles River_1_1"))
C11_E2 <- pcoa_duplicates_time_context(pcoa_data_duplicates %>% filter(id == 'C11D8E2'| id == 'C11D9E1' | mouse_id == "2_Charles River_1_1"))
C11_E1_rmD8E2 <- pcoa_duplicates_time_context(pcoa_data_duplicates %>% filter(id == 'C11D9E1' | mouse_id == "1_Charles River_1_1"))
plot_grid(C11D8E2vD9E1, C11_E1_rmD8E2, C11_E1, C11_E2, 
          labels = c("Duplicates", "Exp.1 minus D8E2 sample", "Over time Exp.1", "Over time Exp.2"),
          ncol = 2, label_x = 0, label_y = 1)+
  ggsave("exploratory/notebook/duplicates/C11.pdf", width = 8.5, height = 11)

#C21 Mouse----
C21_E1 <- pcoa_duplicates_time_context(pcoa_data_duplicates %>% filter(id == 'C21D8E2'| id == 'C21D9E1' | mouse_id == "1_Charles River_2_1"))
C21_E2 <- pcoa_duplicates_time_context(pcoa_data_duplicates %>% filter(id == 'C21D8E2'| id == 'C21D9E1' | mouse_id == "2_Charles River_2_1"))
C21_E1_rmD8E2 <- pcoa_duplicates_time_context(pcoa_data_duplicates %>% filter(id == 'C21D9E1' | mouse_id == "1_Charles River_2_1"))
plot_grid(C11D8E2vD9E1, C21_E1_rmD8E2, C21_E1, C21_E2, 
          labels = c("Duplicates", "Exp.1 minus D8E2 sample", "Over time Exp.1", "Over time Exp.2"),
          ncol = 2, label_x = 0, label_y = 1)+
  ggsave("exploratory/notebook/duplicates/C21.pdf", width = 8.5, height = 11)

#C22 Mouse----
C22_E1 <- pcoa_duplicates_time_context(pcoa_data_duplicates %>% filter(id == 'C22D8E2'| id == 'C22D9E1' | mouse_id == "1_Charles River_2_2"))
C22_E2 <- pcoa_duplicates_time_context(pcoa_data_duplicates %>% filter(id == 'C22D8E2'| id == 'C22D9E1' | mouse_id == "2_Charles River_2_2"))
C22_E1_rmD8E2 <- pcoa_duplicates_time_context(pcoa_data_duplicates %>% filter(id == 'C22D9E1' | mouse_id == "1_Charles River_2_2"))
plot_grid(C22D8E2vD9E1, C22_E1_rmD8E2, C22_E1, C22_E2, 
          labels = c("Duplicates", "Exp.1 minus D8E2 sample", "Over time Exp.1", "Over time Exp.2"),
          ncol = 2, label_x = 0, label_y = 1)+
  ggsave("exploratory/notebook/duplicates/C22.pdf", width = 8.5, height = 11)

#E11D8E2 dropped during rarefaction, but should not exist because it was found dead prior to D3

#E12 Mouse----
E12_E1 <- pcoa_duplicates_time_context(pcoa_data_duplicates %>% filter(id == 'E12D8E2'| id == 'E12D9E1'| mouse_id == "1_Envigo_1_2"))
E12_E2 <- pcoa_duplicates_time_context(pcoa_data_duplicates %>% filter(id == 'E12D8E2'| id == 'E12D9E1'| mouse_id == "2_Envigo_1_2"))
E12_E1_rmD8E2 <- pcoa_duplicates_time_context(pcoa_data_duplicates %>% filter(id == 'E12D9E1' | mouse_id == "1_Envigo_1_2"))
plot_grid(E12D8E2vD9E1, E12_E1_rmD8E2, E12_E1, E12_E2, 
          labels = c("Duplicates", "Exp.1 minus D8E2 sample", "Over time Exp.1", "Over time Exp.2"),
          ncol = 2, label_x = 0, label_y = 1)+
  ggsave("exploratory/notebook/duplicates/E12.pdf", width = 8.5, height = 11)

#E21 Mouse----
E21_E1 <- pcoa_duplicates_time_context(pcoa_data_duplicates %>% filter(id == 'E21D8E2'| id == 'E21D9E1'| mouse_id == "1_Envigo_2_1"))
E21_E2 <- pcoa_duplicates_time_context(pcoa_data_duplicates %>% filter(id == 'E21D8E2'| id == 'E21D9E1'| mouse_id == "2_Envigo_2_1"))
E21_E1_rmD8E2 <- pcoa_duplicates_time_context(pcoa_data_duplicates %>% filter(id == 'E21D9E1' | mouse_id == "1_Envigo_2_1"))
plot_grid(E21D8E2vD9E1, E21_E1_rmD8E2, E21_E1, E21_E2, 
          labels = c("Duplicates", "Exp.1 minus D8E2 sample", "Over time Exp.1", "Over time Exp.2"),
          ncol = 2, label_x = 0, label_y = 1)+
  ggsave("exploratory/notebook/duplicates/E21.pdf", width = 8.5, height = 11)

#E22 Mouse----
E22_E1 <- pcoa_duplicates_time_context(pcoa_data_duplicates %>% filter(id == 'E22D8E2'| id == 'E22D9E1'| mouse_id == "1_Envigo_2_2"))
E22_E2 <- pcoa_duplicates_time_context(pcoa_data_duplicates %>% filter(id == 'E22D8E2'| id == 'E22D9E1'| mouse_id == "2_Envigo_2_2"))
E22_E1_rmD8E2 <- pcoa_duplicates_time_context(pcoa_data_duplicates %>% filter(id == 'E22D9E1' | mouse_id == "1_Envigo_2_2"))
plot_grid(E22D8E2vD9E1, E22_E1_rmD8E2, E22_E1, E22_E2, 
          labels = c("Duplicates", "Exp.1 minus D8E2 sample", "Over time Exp.1", "Over time Exp.2"),
          ncol = 2, label_x = 0, label_y = 1)+
  ggsave("exploratory/notebook/duplicates/E22.pdf", width = 8.5, height = 11)

#J11 Mouse----
J11_E1 <- pcoa_duplicates_time_context(pcoa_data_duplicates %>% filter(id == 'J11D8E2'| id == 'J11D9E1'| mouse_id == "1_Jackson_1_1"))
J11_E2 <- pcoa_duplicates_time_context(pcoa_data_duplicates %>% filter(id == 'J11D8E2'| id == 'J11D9E1'| mouse_id == "2_Jackson_1_1"))
J11_E1_rmD8E2 <- pcoa_duplicates_time_context(pcoa_data_duplicates %>% filter(id == 'J11D9E1' | mouse_id == "1_Jackson_1_1"))
plot_grid(J11D8E2vD9E1, J11_E1_rmD8E2, J11_E1, J11_E2, 
          labels = c("Duplicates", "Exp.1 minus D8E2 sample", "Over time Exp.1", "Over time Exp.2"),
          ncol = 2, label_x = 0, label_y = 1)+
  ggsave("exploratory/notebook/duplicates/J11.pdf", width = 8.5, height = 11)

#J21D8E2 should not exist because it was found dead prior to Day 3
#J21 Mouse----
J21_E1 <- pcoa_duplicates_time_context(pcoa_data_duplicates %>% filter(id == 'J21D8E2'| id == 'J21D9E1'| mouse_id == "1_Jackson_2_1"))
J21_E2 <- pcoa_duplicates_time_context(pcoa_data_duplicates %>% filter(id == 'J21D8E2'| id == 'J21D9E1'| mouse_id == "2_Jackson_2_1"))
J21_E1_rmD8E2 <- pcoa_duplicates_time_context(pcoa_data_duplicates %>% filter(id == 'J21D9E1' | mouse_id == "1_Jackson_2_1"))
plot_grid(J21D8E2vD9E1, J21_E1_rmD8E2, J21_E1, J21_E2, 
          labels = c("Duplicates", "Exp.1 minus D8E2 sample", "Over time Exp.1", "Over time Exp.2"),
          ncol = 2, label_x = 0, label_y = 1)+
  ggsave("exploratory/notebook/duplicates/J21.pdf", width = 8.5, height = 11)

#J22 Mouse----
J22_E1 <- pcoa_duplicates_time_context(pcoa_data_duplicates %>% filter(id == 'J22D8E2' | id == 'J22D9E1'| mouse_id == "1_Jackson_2_2"))
J22_E2 <- pcoa_duplicates_time_context(pcoa_data_duplicates %>% filter(id == 'J22D8E2' | id == 'J22D9E1'| mouse_id == "2_Jackson_2_2"))
J22_E1_rmD8E2 <- pcoa_duplicates_time_context(pcoa_data_duplicates %>% filter(id == 'J22D9E1' | mouse_id == "1_Jackson_2_2"))
plot_grid(J22D8E2vD9E1, J22_E1_rmD8E2, J22_E1, J22_E2, 
          labels = c("Duplicates", "Exp.1 minus D8E2 sample", "Over time Exp.1", "Over time Exp.2"),
          ncol = 2, label_x = 0, label_y = 1)+
  ggsave("exploratory/notebook/duplicates/J22.pdf", width = 8.5, height = 11)

#S12 Mouse----
S12_E1 <- pcoa_duplicates_time_context(pcoa_data_duplicates %>% filter(id == 'S12D8E2' | id == 'S12D9E1'| mouse_id == "1_Schloss_1_2"))
S12_E2 <- pcoa_duplicates_time_context(pcoa_data_duplicates %>% filter(id == 'S12D8E2' | id == 'S12D9E1'| mouse_id == "2_Schloss_1_2"))
S12_E1_rmD8E2 <- pcoa_duplicates_time_context(pcoa_data_duplicates %>% filter(id == 'S12D9E1' | mouse_id == "1_Schloss_1_2"))
plot_grid(S12D8E2vD9E1, S12_E1_rmD8E2, S12_E1, S12_E2, 
          labels = c("Duplicates", "Exp.1 minus D8E2 sample", "Over time Exp.1", "Over time Exp.2"),
          ncol = 2, label_x = 0, label_y = 1)+
  ggsave("exploratory/notebook/duplicates/S12.pdf", width = 8.5, height = 11)

#S21 Mouse----
S21_E1 <- pcoa_duplicates_time_context(pcoa_data_duplicates %>% filter(id == 'S21D8E2' | id == 'S21D9E1'| mouse_id == "1_Schloss_2_1"))
S21_E2 <- pcoa_duplicates_time_context(pcoa_data_duplicates %>% filter(id == 'S21D8E2' | id == 'S21D9E1'| mouse_id == "2_Schloss_2_1"))
S21_E1_rmD8E2 <- pcoa_duplicates_time_context(pcoa_data_duplicates %>% filter(id == 'S21D9E1' | mouse_id == "1_Schloss_2_1"))
plot_grid(S21D8E2vD9E1, S21_E1_rmD8E2, S21_E1, S21_E2, 
          labels = c("Duplicates", "Exp.1 minus D8E2 sample", "Over time Exp.1", "Over time Exp.2"),
          ncol = 2, label_x = 0, label_y = 1)+
  ggsave("exploratory/notebook/duplicates/S21.pdf", width = 8.5, height = 11)

#S22 Mouse----
S22_E1 <- pcoa_duplicates_time_context(pcoa_data_duplicates %>% filter(id == 'S22D8E2' | id == 'S22D9E1'| mouse_id == "1_Schloss_2_2"))
S22_E2 <- pcoa_duplicates_time_context(pcoa_data_duplicates %>% filter(id == 'S22D8E2' | id == 'S22D9E1'| mouse_id == "2_Schloss_2_2"))
S22_E1_rmD8E2 <- pcoa_duplicates_time_context(pcoa_data_duplicates %>% filter(id == 'S22D9E1' | mouse_id == "1_Schloss_2_2"))
plot_grid(S22D8E2vD9E1, S22_E1_rmD8E2, S22_E1, S22_E2, 
          labels = c("Duplicates", "Exp.1 minus D8E2 sample", "Over time Exp.1", "Over time Exp.2"),
          ncol = 2, label_x = 0, label_y = 1)+
  ggsave("exploratory/notebook/duplicates/S22.pdf", width = 8.5, height = 11)

#T11 Mouse----
T11_E1 <- pcoa_duplicates_time_context(pcoa_data_duplicates %>% filter(id == 'T11D8E2' | id == 'T11D9E1'| mouse_id == "1_Taconic_1_1"))
T11_E2 <- pcoa_duplicates_time_context(pcoa_data_duplicates %>% filter(id == 'T11D8E2' | id == 'T11D9E1'| mouse_id == "2_Taconic_1_1"))
T11_E1_rmD8E2 <- pcoa_duplicates_time_context(pcoa_data_duplicates %>% filter(id == 'T11D9E1' | mouse_id == "1_Taconic_1_1"))
plot_grid(T11D8E2vD9E1, T11_E1_rmD8E2, T11_E1, T11_E2, 
          labels = c("Duplicates", "Exp.1 minus D8E2 sample", "Over time Exp.1", "Over time Exp.2"),
          ncol = 2, label_x = 0, label_y = 1)+
  ggsave("exploratory/notebook/duplicates/T11.pdf", width = 8.5, height = 11)

#T12 Mouse----
T12_E1 <- pcoa_duplicates_time_context(pcoa_data_duplicates %>% filter(id == 'T12D8E2' | id == 'T12D9E1'| mouse_id == "1_Taconic_1_2"))
T12_E2 <- pcoa_duplicates_time_context(pcoa_data_duplicates %>% filter(id == 'T12D8E2' | id == 'T12D9E1'| mouse_id == "2_Taconic_1_2"))
T12_E1_rmD8E2 <- pcoa_duplicates_time_context(pcoa_data_duplicates %>% filter(id == 'T12D9E1' | mouse_id == "1_Taconic_1_2"))
plot_grid(T12D8E2vD9E1, T12_E1_rmD8E2, T12_E1, T12_E2, 
          labels = c("Duplicates", "Exp.1 minus D8E2 sample", "Over time Exp.1", "Over time Exp.2"),
          ncol = 2, label_x = 0, label_y = 1)+
  ggsave("exploratory/notebook/duplicates/T12.pdf", width = 8.5, height = 11)

#T21 Mouse----
T21_E1 <- pcoa_duplicates_time_context(pcoa_data_duplicates %>% filter(id == 'T21D8E2' | id == 'T21D9E1'| mouse_id == "1_Taconic_2_1"))
T21_E2 <- pcoa_duplicates_time_context(pcoa_data_duplicates %>% filter(id == 'T21D8E2' | id == 'T21D9E1'| mouse_id == "2_Taconic_2_1"))
T21_E1_rmD8E2 <- pcoa_duplicates_time_context(pcoa_data_duplicates %>% filter(id == 'T21D9E1' | mouse_id == "1_Taconic_2_1"))
plot_grid(T21D8E2vD9E1, T21_E1_rmD8E2, T21_E1, T21_E2, 
          labels = c("Duplicates", "Exp.1 minus D8E2 sample", "Over time Exp.1", "Over time Exp.2"),
          ncol = 2, label_x = 0, label_y = 1)+
  ggsave("exploratory/notebook/duplicates/T21.pdf", width = 8.5, height = 11)

#T22 was NA in Experiment 1 but should have been present in Experiment 2
#Y11 is missing for D8 Experiment 2
#Y12 Mouse----
Y12_E1 <- pcoa_duplicates_time_context(pcoa_data_duplicates %>% filter(id == 'Y12D8E2' | id == 'Y12D9E1'| mouse_id == "1_Young_1_2"))
Y12_E2 <- pcoa_duplicates_time_context(pcoa_data_duplicates %>% filter(id == 'Y12D8E2' | id == 'Y12D9E1'| mouse_id == "2_Young_1_2"))
Y12_E1_rmD8E2 <- pcoa_duplicates_time_context(pcoa_data_duplicates %>% filter(id == 'Y12D9E1' | mouse_id == "1_Young_1_2"))
plot_grid(Y12D8E2vD9E1, Y12_E1_rmD8E2, Y12_E1, Y12_E2, 
          labels = c("Duplicates", "Exp.1 minus D8E2 sample", "Over time Exp.1", "Over time Exp.2"),
          ncol = 2, label_x = 0, label_y = 1)+
  ggsave("exploratory/notebook/duplicates/Y12.pdf", width = 8.5, height = 11)

#Y13 Mouse----
#Y13 should only exist in Experiment 1
Y13_E1 <- pcoa_duplicates_time_context(pcoa_data_duplicates %>% filter(id == 'Y13D8E2' | id == 'Y13D9E1'| mouse_id == "1_Young_1_3"))
Y13_E2 <- pcoa_duplicates_time_context(pcoa_data_duplicates %>% filter(id == 'Y13D8E2' | id == 'Y13D9E1'| mouse_id == "2_Young_1_3"))
Y13_E1_rmD8E2 <- pcoa_duplicates_time_context(pcoa_data_duplicates %>% filter(id == 'Y13D9E1' | mouse_id == "1_Young_1_3"))
plot_grid(Y13D8E2vD9E1, Y13_E1_rmD8E2, Y13_E1, Y13_E2, 
          labels = c("Duplicates", "Exp.1 minus D8E2 sample", "Over time Exp.1", "Over time Exp.2"),
          ncol = 2, label_x = 0, label_y = 1)+
  ggsave("exploratory/notebook/duplicates/Y13.pdf", width = 8.5, height = 11)

#Y21 Mouse----
Y21_E1 <- pcoa_duplicates_time_context(pcoa_data_duplicates %>% filter(id == 'Y21D8E2' | id == 'Y21D9E1'| mouse_id == "1_Young_2_1"))
Y21_E2 <- pcoa_duplicates_time_context(pcoa_data_duplicates %>% filter(id == 'Y21D8E2' | id == 'Y21D9E1'| mouse_id == "2_Young_2_1"))
Y21_E1_rmD8E2 <- pcoa_duplicates_time_context(pcoa_data_duplicates %>% filter(id == 'Y21D9E1' | mouse_id == "1_Young_2_1"))
plot_grid(Y21D8E2vD9E1, Y21_E1_rmD8E2, Y21_E1, Y21_E2, 
          labels = c("Duplicates", "Exp.1 minus D8E2 sample", "Over time Exp.1", "Over time Exp.2"),
          ncol = 2, label_x = 0, label_y = 1)+
  ggsave("exploratory/notebook/duplicates/Y21.pdf", width = 8.5, height = 11)

#Y22 Mouse----
Y22D8E2vD9E1 <- plot_pcoa_duplicates(pcoa_data_duplicates %>% filter(id == 'Y22D8E2' | id == 'Y22D9E1'))
Y22_E1 <- pcoa_duplicates_time_context(pcoa_data_duplicates %>% filter(id == 'Y22D8E2' | id == 'Y22D9E1'| mouse_id == "1_Young_2_2"))
Y22_E2 <- pcoa_duplicates_time_context(pcoa_data_duplicates %>% filter(id == 'Y22D8E2' | id == 'Y22D9E1'| mouse_id == "2_Young_2_2"))
Y22_E1_rmD8E2 <- pcoa_duplicates_time_context(pcoa_data_duplicates %>% filter(id == 'Y22D9E1' | mouse_id == "1_Young_2_2"))
plot_grid(Y22D8E2vD9E1, Y22_E1_rmD8E2, Y22_E1, Y22_E2, 
          labels = c("Duplicates", "Exp.1 minus D8E2 sample", "Over time Exp.1", "Over time Exp.2"),
          ncol = 2, label_x = 0, label_y = 1)+
  ggsave("exploratory/notebook/duplicates/Y22.pdf", width = 8.5, height = 11)
