#Source: https://github.com/SchlossLab/ML_pipeline_microbiome/blob/master/test/code/prepare_input_data.R

#Script prepares the dataset to have C. difficile colonization status on day 7 post-challenge as the outcome. The features are Otu abundances at Day -1, 0, and 1

library(tidyverse)

######################## DATA PREPARATION #############################
# Features: 16S rRNA gene sequences(OTUs) in the stool on Day -1 of the experiment
# Labels: - C. difficile colonization status in the stool 7 days post-challenge

# Read in metadata and select only sample id, day, and cfu_d5 columns
source("code/functions.R")
metadata <- metadata %>% select(id, day, clearance_status_d7)
# Read in OTU table and remove label and numOtus columns
shared <- read.delim('data/process/vendors.subsample.shared', header=T, sep='\t') %>%
  select(-label, -numOtus)
# Function that will do the following
# Merge metadata and OTU table.
# Filter merged dataframe to just timepoints we want to use community structure to predict C. difficile CFU 7 days post-challenge
# Then remove the sample ID and day columns
select_timepoint <- function(timepoint){
  data <- inner_join(metadata, shared, by=c("id"="Group")) %>% 
    filter(day == timepoint) %>% 
    select(-id, -day) %>%
    drop_na() %>%
    select(clearance_status_d7, everything()) %>%
    rename(dx=clearance_status_d7) %>% #rename column, the clearance_status_d7 values are what we want to predict
    select(-Otu0020) %>% #remove C. difficile (Otu0020) from the input data since C. difficile colonization status at day 7 is what we're predicting
    write_csv(paste0("data/process/classification_input_day", timepoint, "_data.csv"))
}
class_data_dn1 <- select_timepoint(-1)
class_data_d1 <- select_timepoint(1)
class_data_d0 <- select_timepoint(0)

###################################################################
