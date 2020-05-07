source("code/functions.R")

#Interpretation----
# Source: https://github.com/SchlossLab/Topcuoglu_ML_XXX_2019/blob/master/code/learning/get_feature_rankings.R
# https://github.com/SchlossLab/Topcuoglu_ML_XXX_2019/blob/master/code/merge_feature_ranks.sh

#Starting from combined_all_imp_features_cor_results file----

# Add column to capture data_split number
data_split = as.numeric(c(0:99))
create_feature_rankings <- function(data){
  # If the models are linear, we saved the weights of every OTU for each datasplit
  # We want to plot the ranking of families for linear models.
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

# Correlated files for linear models has the weights for families in trained model combined from 100 data splits:

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
    # We want to use/plot only the top 20 highest ranked families
    # Initial step is to get which are the highest 20 ranked families by looking at their median rank
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
    
    # Here we want to only grab the data (rank info from 100 datasplits) of only the top 20 median ranked families
    # The imp data will be returned for Figure 3 where we plot each rank info for each data-split of the 20 top families
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
#     2. Returns the top 20 ranked (1-20 with 1 being the highest ranked) families (ranks of the OTU for 100 splits)
get_feature_ranked_files <- function(file_name, model_name){
  importance_data <- read_tsv(file_name)
  ranks <- get_interp_info(importance_data, model_name) %>%
    as.data.frame()
  return(ranks)
}

# This function:
#     1. Top 20 ranked (1-20 lowest rank) families (ranks of the family for 100 splits)
#     2. Returns a plot. Each datapoint is the rank of the family at one datasplit.

plot_feature_ranks <- function(data){
  # Plot from highest median ranked family to least (only top 5) and thir ranks that lay between 1-100
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
#     1. Grabs the taxonomy information for all families
#     2. Grabs the top 20 important taxa and their feature rankings
#     3. Merges the 2 tables (taxa and importance) by the family name
#     4. Make a new column that merges taxa info and family name
#     5. Turns that to a list to use as labels in the plot

get_taxa_info_as_labels <- function(data){
  
  # Grab the names of the top 20 families in the order of their median rank  
  families <- data %>% 
    group_by(key) %>% 
    summarise(imp = median(rank)) %>% 
    arrange(-imp) 
  # Names of the y-axis labels
  taxa_info <- read.delim('data/process/vendors.family.taxonomy', header=T, sep='\t') %>% 
    select(-Size) %>% 
    mutate(key=OTU) %>% 
    select(-OTU)
  
  taxa_families <- inner_join(families, taxa_info, by="key") %>% 
    mutate(taxa=Taxonomy) %>% 
    select(taxa, imp)
  
  return(taxa_families$taxa)
}

######################################################################
#--------------Run the functions and plot feature ranks ----------#----
######################################################################

#Analysis of 3 classification models using 60:40 split----
dn1_60_rank <- create_rank_file("data/process/classification/combined_all_imp_features_cor_results_dayn1_60_family.csv", "dayn1_60_family")
d0_60_rank <- create_rank_file("data/process/classification/combined_all_imp_features_cor_results_day0_60_family.csv", "day0_60_family")
d1_60_rank <- create_rank_file("data/process/classification/combined_all_imp_features_cor_results_day1_60_family.csv", "day1_60_family")

interp_D1_60 <- get_feature_ranked_files("data/process/classification/day1_60_family_L2_Logistic_Regression_feature_ranking.tsv", "L2_Logistic_Regression")
interp_graph_D1_60 <- plot_feature_ranks(interp_D1_60) +
  scale_x_discrete(name = "Day 1 Community",
                   labels = get_taxa_info_as_labels(interp_D1_60)) +
  theme(axis.text.x=element_text(size = 12, colour='black'))
# -------------------------------------------------------------------->


interp_Dn1_60 <- get_feature_ranked_files("data/process/classification/dayn1_60_family_L2_Logistic_Regression_feature_ranking.tsv", "L2_Logistic_Regression")
interp_graph_Dn1_60 <- plot_feature_ranks(interp_Dn1_60) +
  scale_x_discrete(name = "Day -1 Community",
                   labels = get_taxa_info_as_labels(interp_Dn1_60)) +
  theme(axis.text.x=element_text(size = 12, colour='black'))

interp_D0_60 <- get_feature_ranked_files("data/process/classification/day0_60_family_L2_Logistic_Regression_feature_ranking.tsv", "L2_Logistic_Regression")
interp_graph_D0_60 <- plot_feature_ranks(interp_D0_60) +
  scale_x_discrete(name = "Day 0 Community",
                   labels = get_taxa_info_as_labels(interp_D0_60)) +
  theme(axis.text.x=element_text(size = 12, colour='black'))

######################################################################
#-----------------------Save figure as .pdf ------------------------ #
######################################################################
linear <- plot_grid(interp_graph_D0_60, interp_graph_Dn1_60, interp_graph_D1_60, labels = c("A", "B", "C"), align = 'v', ncol = 1)

ggdraw(add_sub(linear, "Feature ranks", vpadding=grid::unit(0,"lines"), y=5, x=0.7, vjust=4.75, size=15))

ggsave("results/figures/class_interp_60_family.pdf", plot = linear, width = 6, height = 9.2, dpi=300)


#Function to find which significant families are shared 
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
# 3 Families: "Bacteroidaceae"    "Coriobacteriaceae" "Ruminococcaceae"  
#Shared between -1 & 0:
dn1_0_shared <- intersect_all(dn1_taxa, d0_taxa)
# 5 families: [1] "Bacteroidaceae"     "Coriobacteriaceae" 
# "Ruminococcaceae"    "Staphylococcaceae" 
# "Deferribacteraceae" 
#Shared between 0 & 1:
d1_0_shared <- intersect_all(d0_taxa, d1_taxa)
#7 families: 
# "Lachnospiraceae"              "Coriobacteriaceae"           
# "Verrucomicrobiaceae"          "Enterococcaceae"             
# "Bifidobacteriaceae"           "Lactobacillales_unclassified"
# "Ruminococcaceae"              "Bacteroidaceae"
#Shared between -1 & 1:
dn1_1_shared <- intersect_all(dn1_taxa, d1_taxa)
#5 families: "Bacteroidaceae", "Coriobacteriaceae", "Peptostreptococcaceae", "Ruminococcaceae", "Proteobacteria_unclassified" "Clostridiaceae_1" 

#Write list of important families from the 3 models to .csvs to use in plotting taxa
tibble(dayn1_interp_families = dn1_taxa ) %>% 
  write_csv(paste0("data/process/interp_families_dn1.csv"))
tibble(day0_interp_families = d0_taxa ) %>% 
  write_csv(paste0("data/process/interp_families_d0.csv"))
tibble(day1_interp_families = d1_taxa ) %>% 
  write_csv(paste0("data/process/interp_families_d1.csv"))
