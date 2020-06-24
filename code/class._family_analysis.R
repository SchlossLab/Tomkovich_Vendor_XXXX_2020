source("code/functions.R")

#Adapted from:
#https://github.com/SchlossLab/Topcuoglu_ML_XXX_2019/blob/master/code/learning/Figure2.R
# https://github.com/SchlossLab/Sze_SCFACRC_mBio_2019/blob/master/code/plot_classification_regression.R
# https://github.com/SchlossLab/Sze_SCFACRC_mBio_2019/blob/master/code/plot_classification_fit.R
# https://github.com/SchlossLab/Sze_SCFACRC_mBio_2019/blob/master/code/plot_regression_fit.R


#Classification: L2 Logistic Regression model used to predict C. difficile colonization status on Day 7 post-challenge (cleared/not_detectable), colonized)----
#Based on the family relative abundances on Day -1, 0 or 1 of the experiment (C. difficile OTU 20 was removed prior to pooling OTU relative abundance data into families)

#Trained and tested with 60/40 splits----
#Read in the cvAUCs and test AUCs for 100 splits, need to add a column indicating input Otu data that was used
dayn1_file <- list.files(path= 'data/process/classification', pattern='combined_best_hp_results_dayn1_60_family', full.names = TRUE) 
dayn1 <- map_df(dayn1_file, read_files) %>% 
  mutate(input_data = "dayn1") 
day0_file <- list.files(path= 'data/process/classification', pattern='combined_best_hp_results_day0_60_family', full.names = TRUE) 
day0 <- map_df(day0_file, read_files) %>% 
  mutate(input_data = "day0") 
day1_file <- list.files(path= 'data/process/classification', pattern='combined_best_hp_results_day1_60_family', full.names = TRUE) 
day1 <- map_df(day1_file, read_files) %>% 
  mutate(input_data = "day1") 
#Combine the 3 dataframes 
best_performance <- rbind(dayn1, day0, day1) %>% 
  melt_data() 
#Median AUC (based on cvAUCs and testAUCs)  
median_performance <- best_performance %>% 
  group_by(input_data) %>% summarise(median = median(AUC), n = n())

#Plot the cvAUC and testAUC of 100 datasplits----
performance <- ggplot(best_performance, aes(x = fct_reorder(input_data, AUC), y = AUC, fill = Performance)) +
  geom_boxplot(alpha=0.5, fatten = 4) +
  geom_hline(yintercept = 0.5, linetype="dashed") +
  scale_fill_manual(values=c("blue4", "springgreen4"), 
                    guide = guide_legend(reverse = TRUE)) +
  coord_flip() +
  scale_y_continuous(name = "AUROC",
                     breaks = seq(0.4, 1, 0.1),
                     limits=c(0.4, 1),
                     expand=c(0,0)) +
 scale_x_discrete(name = "",
                  labels=c("Day 1",
                           "Day -1",
                           "Day 0")) +
  theme_bw() +
  ggtitle("Family level")+ #Title plot
  theme(plot.title = element_text(hjust = 0.5)) #Center plot titile
  # theme(plot.margin=unit(c(1.1,1.1,1.1,1.1),"cm"),
  #       legend.justification=c(1,0),
  #       legend.position=c(1,0),
  #       #legend.position="bottom",
  #       legend.title = element_blank(),
  #       legend.background = element_rect(linetype="solid", color="black", size=0.5),
  #       legend.box.margin=margin(c(12,12,12, 12)),
  #       legend.text=element_text(size=18),
  #       #legend.title=element_text(size=22),
  #       panel.grid.major.y = element_blank(),
  #       panel.grid.major.x = element_line( size=0.6),
  #       panel.grid.minor = element_blank(),
  #       panel.background = element_blank(),
  #       text = element_text(size = 12),
  #       axis.text.x=element_text(size = 20, colour='black'),
  #       axis.text.y=element_text(size = 20, colour='black'),
  #       axis.title.y=element_text(size = 24),
  #       axis.title.x=element_text(size = 24),
  #       panel.border = element_rect(linetype="solid", colour = "black", fill=NA, size=1.5))


######################################################################
#-----------------------Save figure as .png ------------------------ #
######################################################################

save_plot("results/figures/class._logistic_regression_60_family.png", plot = performance, base_aspect_ratio = 1.8)

