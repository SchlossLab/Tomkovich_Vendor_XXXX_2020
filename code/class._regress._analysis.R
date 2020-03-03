source("code/functions.R")

#Adapted from:
#https://github.com/SchlossLab/Topcuoglu_ML_XXX_2019/blob/master/code/learning/Figure2.R
# https://github.com/SchlossLab/Sze_SCFACRC_mBio_2019/blob/master/code/plot_classification_regression.R
# https://github.com/SchlossLab/Sze_SCFACRC_mBio_2019/blob/master/code/plot_classification_fit.R
# https://github.com/SchlossLab/Sze_SCFACRC_mBio_2019/blob/master/code/plot_regression_fit.R


#Classification: L2 Logistic Regression model used to predict C. difficile colonization status on Day 7 post-challenge (cleared, colonized)----
#Trained with data from 32 mice: 15 cleared, 17 colonized on Day 7
#Based on the Otu relative abundances on Day -1, 0 or 1 of the experiment

#Read in the cvAUCs and test AUCs for 100 splits, need to add a column indicating input Otu data that was used
dayn1_file <- list.files(path= 'data/process/classification', pattern='combined_best_hp_results_dayn1', full.names = TRUE) 
dayn1 <- map_df(dayn1_file, read_files) %>% 
  mutate(input_data = "dayn1") 
day0_file <- list.files(path= 'data/process/classification', pattern='combined_best_hp_results_day0', full.names = TRUE) 
day0 <- map_df(day0_file, read_files) %>% 
  mutate(input_data = "day0") 
day1_file <- list.files(path= 'data/process/classification', pattern='combined_best_hp_results_day1', full.names = TRUE) 
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
  theme_bw() #+
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
#-----------------------Save figure as .pdf ------------------------ #
######################################################################

ggsave("results/figures/class._logistic_regression.pdf", plot = performance)


#Regression: L2 regularized L2-loss support vector regression used to predict C. difficile colonization levels on Day 5 post-challenge----
#Trained with data from 36 mice
#Based on the Otu relative abundances on Day -1, 0 or 1 of the experiment

#Read in the cv(train) and test Rsquared values for 100 splits, need to add a column indicating input Otu data that was used
r_dayn1_file <- list.files(path= 'data/process/regression', pattern='combined_optimum_cost__dayn1.csv', full.names = TRUE) 
r_dayn1 <- map_df(r_dayn1_file, read_files) %>% 
  mutate(input_data = "dayn1") 
r_day0_file <- list.files(path= 'data/process/regression', pattern='combined_optimum_cost__day0.csv', full.names = TRUE) 
r_day0 <- map_df(r_day0_file, read_files) %>% 
  mutate(input_data = "day0") 
r_day1_file <- list.files(path= 'data/process/regression', pattern='combined_optimum_cost__day1.csv', full.names = TRUE) 
r_day1 <- map_df(r_day1_file, read_files) %>% 
  mutate(input_data = "day1") 
#Combine the 3 dataframes 
regress_best_performance <- rbind(r_dayn1, r_day0, r_day1) %>% 
  # Instead of 2 columns with names cv/train_rsquared and test_rsquared
  # We will have 1 column with name Performance that tells us if test or cv(train)
  melt(measure.vars=c('train_Rsquared', 'test_Rsquared')) %>%
  rename(Rsquared=value) %>%
  mutate(Performance = case_when(variable == "train_Rsquared" ~ 'Cross-validation', variable == "test_Rsquared" ~ 'Testing')) %>%
  group_by(Performance) 
#Median Rsquared (based on Rsquared and Rsquared)  
median_regress_performance <- regress_best_performance %>% 
  group_by(input_data) %>% summarise(median = median(Rsquared), n = n())

#Plot the cvRsquared and testRsquared of 100 datasplits----
regress_performance <- ggplot(regress_best_performance, aes(x = fct_reorder(input_data, Rsquared), y = Rsquared, fill = Performance)) +
  geom_boxplot(alpha=0.5, fatten = 4) +
  scale_fill_manual(values=c("blue4", "springgreen4"), 
                    guide = guide_legend(reverse = TRUE)) +
  coord_flip() +
  scale_y_continuous(name = "Rsquared",
                     breaks = seq(0, 1, 0.1),
                     limits=c(0, 1),
                     expand=c(0,0)) +
  scale_x_discrete(name = "",
                   labels=c("Day 1",
                            "Day 0",
                            "Day -1")) +
  theme_bw()

######################################################################
#-----------------------Save figure as .pdf ------------------------ #
######################################################################

ggsave("results/figures/regress_l2_loss.pdf", plot = regress_performance)
