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
    theme_classic()

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




