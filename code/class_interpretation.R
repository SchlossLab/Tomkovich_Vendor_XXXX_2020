source("code/functions.R")

#Interpretation----
# Source: https://github.com/SchlossLab/Topcuoglu_ML_XXX_2019/blob/master/code/learning/get_feature_rankings.R
# https://github.com/SchlossLab/Topcuoglu_ML_XXX_2019/blob/master/code/merge_feature_ranks.sh


#After running class_interpretation.R and merge_feature_ranks.sh in ML_pipeline_microbiome

# ------------------- Re-organize feature importance  ----------------->
# This function:
#     1. Takes in a dataframe (different data for each model) and the model name
#     2. If the models are linear, returns the median rank of the top ranked 5 features
#     3. If the models are not linear, returns the permutation importance results for:
#         - Correlated and non-correlated OTUs:
#         - Top 20 features or feature groups will be listed
#         - New AUROC which will differ from original AUROC after permutation
get_interp_info <- function(data, model_name){
  if("key" %in% colnames(data)){
    # If the models are linear, we used get_feature_rankings.R and then merge_feature_ranks.sh first
    # The created file after those 2 steps will be used in this function,
    # Data format is:
    #         The OTU names are in 1 column(repeated for 100 datasplits)
    #         The ranks based on absolute weights are in 1 column(for each of the datasplits)
    #		  The weight value is on another column
    #     The sign of the weight is on another column
    # We want to use/plot only the top 20 highest ranked OTUs
    # Initial step is to get which are the highest 20 ranked OTUs by looking at their median rank
    # 1. We group by OTU name to make sure we are taking all the data-splits into account
    imp_first_20 <- data %>%
      # 2. Group by the OTU name and compute median rank for each OTU
      group_by(key) %>%
      summarise(median_rank = median(rank)) %>%
      # 3. Arrange from highest ranked 1, descending
      arrange(median_rank) %>%
      # 4. Grab only the highest ranked 20
      head(n=20) %>%
      select(key, median_rank)
    
    # Here we want to only grab the data (rank info from 100 datasplits) of only the top 20 median ranked OTUs
    # The imp data will be returned for Figure where we plot each rank info for each data-split of the 20 top OTUs
    imp <- data %>%
      filter(key %in% imp_first_20$key) %>%
      group_by(key)
    
    
    
  }
  # If we want to calculate the permutation importance results for interpretation
  # Then we use the files without the weight information but with permutation results
  else{
    if("names" %in% colnames(data)){ # If the file has non-correlated OTUs
      non_correlated_data <- data %>%
        # 1. Group by the OTU names and calculate median and sd for auc change
        group_by(names) %>%
        summarise(imp = median(new_auc), sd_imp = sd(new_auc))
      # Order the dataframe from smallest new_auc to largest.
      # Because the smallest new_auc means that that OTU decreased AUC a lot when permuted
      imp <- non_correlated_data %>%
        arrange(imp)
    }
    else if("X1" %in% colnames(data)){
      # The file doesn't have "names" column which means these are correlated OTU groups
      # The file has correlated OTUs and their total percent auc change per group in one row
      # Each row has different groups of OTUs that are correlated together
      #     1. We will group by the first OTU (since it is only present in one group only)
      #         This will group all the datasplits for that OTU group together
      #     2. We then get the median percent auc change of that correlated OTU group
      correlated_data <- data %>%
        group_by(X1) %>%
        summarise(imp = median(new_auc), sd_imp = sd(new_auc))
      #     3. We will now only take the first 5 and add the other OTUs to the row.
      #       We have the new_auc for each correlated group of OTUs in a row
      #       We will also have all the OTU names in the group in the same row.
      imp <- correlated_data %>%
        arrange(imp) %>%
        head(5) %>%
        inner_join(data, by="X1") %>% # Add all the other OTUs in the group back to the data
        select(-new_auc, -model)
    }
    else{
      print("linear model")
      imp <- NULL
    }
    
  }
  return(imp)
}
# -------------------------------------------------------------------->

######################################################################
#----------------- Define the functions we will use -----------------#
######################################################################

# ------------------- Re-organize feature importance  ----------------->
# This function:
#     1. Takes in a combined (100 split) feature rankings for each model) and the model name
#     2. Returns the top 20 ranked (1-20 with 1 being the highest ranked) OTUs (ranks of the OTU for 100 splits)
get_feature_ranked_files <- function(file_name, model_name){
  importance_data <- read_tsv(file_name)
  ranks <- get_interp_info(importance_data, model_name) %>%
    as.data.frame()
  return(ranks)
}

# This function:
#     1. Top 20 ranked (1-20 lowest rank) OTUs (ranks of the OTU for 100 splits)
#     2. Returns a plot. Each datapoint is the rank of the OTU at one datasplit.

plot_feature_ranks <- function(data){
  # Plot from highest median ranked OTU to least (only top 20) and thir ranks that lay between 1-100
  # Rank 1 is the highest rank
  plot <- ggplot(data, aes(reorder(data$key, -data$rank, FUN = median), data$rank)) +
    geom_point(aes(colour= factor(data$sign)), size=1.5) + # datapoints lighter color
    scale_color_manual(values=c("#56B4E9","red3", "#999999")) + #blue, red3 and black
    stat_summary(fun.y = function(x) median(x), colour = 'black', geom = "point", size = 3) + # Median darker
    coord_flip(ylim=c(0,100)) +
    theme_classic() +
    theme(plot.margin=unit(c(1.5,3,1.5,3),"mm"),
          legend.position="none",
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.title.y=element_text(size = 11, colour='black', face="bold"), 
          axis.text.y=element_text(size = 8, colour='black'))
  return(plot)
}

# This function:
#     1. Grabs the taxonomy information for all OTUs
#     2. Grabs the top 5 important taxa and their feature rankings
#     3. Merges the 2 tables (taxa and importance) by the OTU name
#     4. Make a new column that merges taxa info and OTU name
#     5. Turns that to a list to use as labels in the plot

get_taxa_info_as_labels <- function(data){
  
  # Grab the names of the top 20 OTUs in the order of their median rank  
  otus <- data %>% 
    group_by(key) %>% 
    summarise(imp = median(rank)) %>% 
    arrange(-imp) 
  # Names of the y-axis labels
  taxa_info <- read.delim('data/process/vendors.taxonomy', header=T, sep='\t') %>% 
    select(-Size) %>% 
    mutate(key=OTU) %>% 
    select(-OTU)
  
  taxa_otus <- inner_join(otus, taxa_info, by="key") %>% 
    mutate_if(is.character, str_to_upper) %>%
    mutate(taxa=gsub("(.*);.*","\\1",Taxonomy)) %>% 
    mutate(taxa=gsub("(.*)_.*","\\1",Taxonomy)) %>% 
    mutate(taxa=gsub("(.*);.*","\\1",Taxonomy)) %>% 
    mutate(taxa=gsub(".*;","",taxa)) %>% 
    mutate(taxa=gsub("(.*)_.*","\\1",taxa)) %>% 
    mutate(taxa=gsub('[0-9]+', '', taxa)) %>% 
    mutate(taxa=str_remove_all(taxa, "[(100)]")) %>% 
    select(key, taxa, imp) %>% 
    unite(key, taxa, key, sep=" (") %>% 
    mutate(key = paste(key,")", sep="")) %>% 
    mutate(key=paste0(gsub('TU0*', 'TU ', key))) 
  
  
  return(taxa_otus$key)
}

######################################################################
#--------------Run the functions and plot feature ranks ----------#----
######################################################################

logit_imp_D1 <- get_feature_ranked_files("data/process/classification/day1/combined_L2_Logistic_Regression_feature_ranking.tsv", "L2_Logistic_Regression")
logit_graph_D1 <- plot_feature_ranks(logit_imp_D1) +
  scale_x_discrete(name = expression(paste(L[2], "-regularized logistic regression")),
                   labels = get_taxa_info_as_labels(logit_imp_D1)) +
  theme(axis.text.x=element_text(size = 12, colour='black'))
# -------------------------------------------------------------------->

# Alternative approach starting from combined_all_imp_features_cor_results file----

# Add column to capture data_split number
data_split = as.numeric(c(0:99))
create_feature_rankings <- function(data){
  # If the models are linear, we saved the weights of every OTU for each datasplit
  # We want to plot the ranking of OTUs for linear models.
  # 1. Get dataframe transformed into long form
  #         The OTU names are in 1 column
  #         The weight value are in 1 column
  weights <- data %>%
    select(-Bias, -model, -data_split) %>%
    gather(factor_key=TRUE) %>%
    mutate(sign = case_when(value<0 ~ "negative",
                            value>0 ~ "positive",
                            value==0 ~ "zero"))
  
  # 2. We change all the weights to their absolute value
  #       Because we want to see which weights are the largest
  weights$value <- abs(weights$value)
  # 3.  a) Order the dataframe from largest weights to smallest.
  #     b) Rank them based on that hierarchy
  #     c) Put the signs back to weights
  #     d) select the OTU names, weights and ranks, return dataframe
  ranks <- weights %>%
    arrange(desc(value)) %>%
    mutate(rank = 1:nrow(weights)) %>%
    mutate(value = case_when(sign=="negative" ~ value*-1,
                             sign=="positive"~ value,
                             sign=="zero" ~ value)) %>%
    select(key, value, rank, sign)
  
  return(ranks)
}


######################################################################
#--------------Run the functions and create ranking files ----------#
######################################################################

# ----------- Read in combined saved weights for linear models in data/process/classification folder ----------> 

# Correlated files for linear models has the weights for OTUs in trained model combined from 100 data splits:

# Take each combined file, add a column specifying data split, group_by data splits, split dataframe into list of 100 data frames, and add ranks as a column to it using lapply
# Then save that file with the "_feature_ranking_#" extension
#Arguments: 
#file_name = combined_all_imp_features_cor_results_$model.csv outputted from pipeline
#day_of_input_data = The experiment day relative abundances that were used to train the model
create_rank_file <- function(file_name, day_of_input_data){ 
  importance_data <- read_files(file_name)
  model_name <- as.character(importance_data$model[1])# get the model name from table
  split_list <- importance_data %>% 
    mutate(data_split = data_split) %>% 
    group_by(data_split) %>% 
    group_split
  data <- lapply(split_list, create_feature_rankings) 
  merged_data <- bind_rows(data[1:100]) 
  write_tsv(merged_data, paste0("data/process/classification/", day_of_input_data, "_L2_Logistic_Regression_feature_ranking.tsv"))
}

dn1_rank <- create_rank_file("data/process/classification/combined_all_imp_features_cor_results_dayn1.csv", "dayn1")
d0_rank <- create_rank_file("data/process/classification/combined_all_imp_features_cor_results_day0.csv", "day0")
d1_rank <- create_rank_file("data/process/classification/combined_all_imp_features_cor_results_day1.csv", "day1")

# ------------------- Re-organize feature importance  ----------------->
# This function:
#     1. Takes in a dataframe (different data for each model) and the model name
#     2. If the models are linear, returns the median rank of the top ranked 20 features
get_interp_info <- function(data, model_name){
    # The created file after those 2 steps will be used in this function,
    # Data format is:
    #         The OTU names are in 1 column(repeated for 100 datasplits)
    #         The ranks based on absolute weights are in 1 column(for each of the datasplits)
    #		  The weight value is on another column
    #     The sign of the weight is on another column
    # We want to use/plot only the top 20 highest ranked OTUs
    # Initial step is to get which are the highest 20 ranked OTUs by looking at their median rank
    # 1. We group by OTU name to make sure we are taking all the data-splits into account
    imp_first_20 <- data %>%
      # 2. Group by the OTU name and compute median rank for each OTU
      group_by(key) %>%
      summarise(median_rank = median(rank)) %>%
      # 3. Arrange from highest ranked 1, descending
      arrange(median_rank) %>%
      # 4. Grab only the highest ranked 20
      head(n=20) %>%
      select(key, median_rank)
    
    # Here we want to only grab the data (rank info from 100 datasplits) of only the top 20 median ranked OTUs
    # The imp data will be returned for Figure 3 where we plot each rank info for each data-split of the 20 top OTUs
    imp <- data %>%
      filter(key %in% imp_first_20$key) %>%
      group_by(key)
    return(imp)
}
  
# -------------------------------------------------------------------->

######################################################################
#----------------- Define the functions we will use -----------------#
######################################################################

# ------------------- Re-organize feature importance  ----------------->
# This function:
#     1. Takes in a combined (100 split) feature rankings for each model) and the model name
#     2. Returns the top 20 ranked (1-20 with 1 being the highest ranked) OTUs (ranks of the OTU for 100 splits)
get_feature_ranked_files <- function(file_name, model_name){
  importance_data <- read_tsv(file_name)
  ranks <- get_interp_info(importance_data, model_name) %>%
    as.data.frame()
  return(ranks)
}

# This function:
#     1. Top 20 ranked (1-20 lowest rank) OTUs (ranks of the OTU for 100 splits)
#     2. Returns a plot. Each datapoint is the rank of the OTU at one datasplit.

plot_feature_ranks <- function(data){
  # Plot from highest median ranked OTU to least (only top 5) and thir ranks that lay between 1-100
  # Rank 1 is the highest rank
  plot <- ggplot(data, aes(reorder(data$key, -data$rank, FUN = median), data$rank)) +
    geom_point(aes(colour= factor(data$sign)), size=1.5) + # datapoints lighter color
    scale_color_manual(values=c("#56B4E9","red3", "#999999")) +
    stat_summary(fun.y = function(x) median(x), colour = 'black', geom = "point", size = 3) + # Median darker
    coord_flip(ylim=c(0,100)) +
    theme_classic() +
    theme(plot.margin=unit(c(1.5,3,1.5,3),"mm"),
          legend.position="none",
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.title.y=element_text(size = 11, colour='black', face="bold"), 
          axis.text.y=element_text(size = 8, colour='black'))+
    theme(axis.text.y = element_text(face = "italic")) #Have the taxa show up as italics
  return(plot)
}

# This function:
#     1. Grabs the taxonomy information for all OTUs
#     2. Grabs the top 20 important taxa and their feature rankings
#     3. Merges the 2 tables (taxa and importance) by the OTU name
#     4. Make a new column that merges taxa info and OTU name
#     5. Turns that to a list to use as labels in the plot

get_taxa_info_as_labels <- function(data){
  
  # Grab the names of the top 20 OTUs in the order of their median rank  
  otus <- data %>% 
    group_by(key) %>% 
    summarise(imp = median(rank)) %>% 
    arrange(-imp) 
  # Names of the y-axis labels
  taxa_info <- read.delim('data/process/vendors.taxonomy', header=T, sep='\t') %>% 
    select(-Size) %>% 
    mutate(key=OTU) %>% 
    select(-OTU)
  
  taxa_otus <- inner_join(otus, taxa_info, by="key") %>% 
    mutate_if(is.character, str_to_upper) %>%
    mutate(taxa=gsub("(.*);.*","\\1",Taxonomy)) %>% 
    mutate(taxa=gsub("(.*)_.*","\\1",Taxonomy)) %>% 
    mutate(taxa=gsub("(.*);.*","\\1",Taxonomy)) %>% 
    mutate(taxa=gsub(".*;","",taxa)) %>% 
    mutate(taxa=gsub("(.*)_.*","\\1",taxa)) %>% 
    mutate(taxa=gsub('[0-9]+', '', taxa)) %>% 
    mutate(taxa=str_remove_all(taxa, "[(100)]")) %>% 
    select(key, taxa, imp) %>% 
    unite(key, taxa, key, sep=" (") %>% 
    mutate(key = paste(key,")", sep="")) %>% 
    mutate(key=paste0(gsub('TU0*', 'TU ', key))) 
  
  
  return(taxa_otus$key)
}

######################################################################
#--------------Run the functions and plot feature ranks ----------#----
######################################################################

interp_D1 <- get_feature_ranked_files("data/process/classification/day1_L2_Logistic_Regression_feature_ranking.tsv", "L2_Logistic_Regression")
interp_graph_D1 <- plot_feature_ranks(interp_D1) +
  scale_x_discrete(name = "Day 1 Community",
                   labels = get_taxa_info_as_labels(interp_D1)) +
  theme(axis.text.x=element_text(size = 12, colour='black'))
# Graph appears equivalent to logit_graph_D1 so will proceed with alternative method (already wrote over temp files for -1, and 0 communities in ML_pipeline_microbiome)
# -------------------------------------------------------------------->


interp_Dn1 <- get_feature_ranked_files("data/process/classification/dayn1_L2_Logistic_Regression_feature_ranking.tsv", "L2_Logistic_Regression")
interp_graph_Dn1 <- plot_feature_ranks(interp_Dn1) +
  scale_x_discrete(name = "Day -1 Community",
                   labels = get_taxa_info_as_labels(interp_Dn1)) +
  theme(axis.text.x=element_text(size = 12, colour='black'))

interp_D0 <- get_feature_ranked_files("data/process/classification/day0_L2_Logistic_Regression_feature_ranking.tsv", "L2_Logistic_Regression")
interp_graph_D0 <- plot_feature_ranks(interp_D0) +
  scale_x_discrete(name = "Day 0 Community",
                   labels = get_taxa_info_as_labels(interp_D0)) +
  theme(axis.text.x=element_text(size = 12, colour='black'))

######################################################################
#-----------------------Save figure as .pdf ------------------------ #
######################################################################
linear <- plot_grid(interp_graph_D0, interp_graph_Dn1, interp_graph_D1, labels = c("A", "B", "C"), align = 'v', ncol = 1)

ggdraw(add_sub(linear, "Feature ranks", vpadding=grid::unit(0,"lines"), y=5, x=0.7, vjust=4.75, size=15))

ggsave("results/figures/class_interp.pdf", plot = linear, width = 6, height = 9.2, dpi=300)

#Analysis of 3 classification models using 50:50 split instead of 80:20----
dn1_50_rank <- create_rank_file("data/process/classification/combined_all_imp_features_cor_results_dayn1_50.csv", "dayn1_50")
d0_50_rank <- create_rank_file("data/process/classification/combined_all_imp_features_cor_results_day0_50.csv", "day0_50")
d1_50_rank <- create_rank_file("data/process/classification/combined_all_imp_features_cor_results_day1_50.csv", "day1_50")

interp_D1_50 <- get_feature_ranked_files("data/process/classification/day1_50_L2_Logistic_Regression_feature_ranking.tsv", "L2_Logistic_Regression")
interp_graph_D1_50 <- plot_feature_ranks(interp_D1_50) +
  scale_x_discrete(name = "Day 1 Community",
                   labels = get_taxa_info_as_labels(interp_D1_50)) +
  theme(axis.text.x=element_text(size = 12, colour='black'))
# -------------------------------------------------------------------->


interp_Dn1_50 <- get_feature_ranked_files("data/process/classification/dayn1_50_L2_Logistic_Regression_feature_ranking.tsv", "L2_Logistic_Regression")
interp_graph_Dn1_50 <- plot_feature_ranks(interp_Dn1_50) +
  scale_x_discrete(name = "Day -1 Community",
                   labels = get_taxa_info_as_labels(interp_Dn1_50)) +
  theme(axis.text.x=element_text(size = 12, colour='black'))

interp_D0_50 <- get_feature_ranked_files("data/process/classification/day0_50_L2_Logistic_Regression_feature_ranking.tsv", "L2_Logistic_Regression")
interp_graph_D0_50 <- plot_feature_ranks(interp_D0_50) +
  scale_x_discrete(name = "Day 0 Community",
                   labels = get_taxa_info_as_labels(interp_D0_50)) +
  theme(axis.text.x=element_text(size = 12, colour='black'))

######################################################################
#-----------------------Save figure as .pdf ------------------------ #
######################################################################
linear <- plot_grid(interp_graph_D0_50, interp_graph_Dn1_50, interp_graph_D1_50, labels = c("A", "B", "C"), align = 'v', ncol = 1)

ggdraw(add_sub(linear, "Feature ranks", vpadding=grid::unit(0,"lines"), y=5, x=0.7, vjust=4.75, size=15))

ggsave("results/figures/class_interp_50.pdf", plot = linear, width = 6, height = 9.2, dpi=300)

#Analysis of 3 classification models using 60:40 split----
dn1_60_rank <- create_rank_file("data/process/classification/combined_all_imp_features_cor_results_dayn1_60.csv", "dayn1_60")
d0_60_rank <- create_rank_file("data/process/classification/combined_all_imp_features_cor_results_day0_60.csv", "day0_60")
d1_60_rank <- create_rank_file("data/process/classification/combined_all_imp_features_cor_results_day1_60.csv", "day1_60")

interp_D1_60 <- get_feature_ranked_files("data/process/classification/day1_60_L2_Logistic_Regression_feature_ranking.tsv", "L2_Logistic_Regression")
interp_graph_D1_60 <- plot_feature_ranks(interp_D1_60) +
  scale_x_discrete(name = "Day 1 Community",
                   labels = get_taxa_info_as_labels(interp_D1_60)) +
  theme(axis.text.x=element_text(size = 12, colour='black'))
# -------------------------------------------------------------------->


interp_Dn1_60 <- get_feature_ranked_files("data/process/classification/dayn1_60_L2_Logistic_Regression_feature_ranking.tsv", "L2_Logistic_Regression")
interp_graph_Dn1_60 <- plot_feature_ranks(interp_Dn1_60) +
  scale_x_discrete(name = "Day -1 Community",
                   labels = get_taxa_info_as_labels(interp_Dn1_60)) +
  theme(axis.text.x=element_text(size = 12, colour='black'))

interp_D0_60 <- get_feature_ranked_files("data/process/classification/day0_60_L2_Logistic_Regression_feature_ranking.tsv", "L2_Logistic_Regression")
interp_graph_D0_60 <- plot_feature_ranks(interp_D0_60) +
  scale_x_discrete(name = "Day 0 Community",
                   labels = get_taxa_info_as_labels(interp_D0_60)) +
  theme(axis.text.x=element_text(size = 12, colour='black'))

######################################################################
#-----------------------Save figure as .pdf ------------------------ #
######################################################################
linear <- plot_grid(interp_graph_D0_60, interp_graph_Dn1_60, interp_graph_D1_60, labels = c("A", "B", "C"), align = 'v', ncol = 1)

ggdraw(add_sub(linear, "Feature ranks", vpadding=grid::unit(0,"lines"), y=5, x=0.7, vjust=4.75, size=15))

ggsave("results/figures/class_interp_60.pdf", plot = linear, width = 6, height = 9.2, dpi=300)


#Function to find which significant otus are shared 
intersect_all <- function(a,b,...){
  Reduce(intersect, list(a,b,...))
}
# Use results from 60:40 splits since that had the best AUROC results
#List of day -1 imp. taxa:
dn1_taxa <- get_taxa_info_as_labels(interp_Dn1_60)
#List of day 0 imp. taxa:
d0_taxa <- get_taxa_info_as_labels(interp_D0_60)
#List of day 1 imp. taxa:
d1_taxa <- get_taxa_info_as_labels(interp_D1_60)

#Shared between all 3 models:
all_shared <- intersect_all(dn1_taxa, d0_taxa, d1_taxa)
# 0 Taxa
#Shared between -1 & 0:
dn1_0_shared <- intersect_all(dn1_taxa, d0_taxa)
#"Ruminococcaceae (OTU 520)"  "Enterobacteriaceae (OTU 1)"
#Shared between 0 & 1:
d1_0_shared <- intersect_all(d0_taxa, d1_taxa)
# "Erysipelotrichaceae (OTU 234)" "Porphyromonadaceae (OTU 7)"   
# "Lactobacillus (OTU 18)"        "Bacteroides (OTU 2)"
#Shared between -1 & 1:
dn1_1_shared <- intersect_all(dn1_taxa, d1_taxa)
#"Enterococcus (OTU 23)"

#Write list of important OTUs from the 3 models to .csvs to use in plotting taxa
tibble(dayn1_interp_otus = dn1_taxa ) %>% 
  write_csv(paste0("data/process/interp_otus_dn1.csv"))
tibble(day0_interp_otus = d0_taxa ) %>% 
  write_csv(paste0("data/process/interp_otus_d0.csv"))
tibble(day1_interp_otus = d1_taxa ) %>% 
  write_csv(paste0("data/process/interp_otus_d1.csv"))

