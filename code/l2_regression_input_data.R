#Begum's suggestion for formatting input data for ML_pipeline_microbiome
#Source: https://github.com/SchlossLab/ML_pipeline_microbiome/blob/master/test/code/prepare_input_data.R

#Script prepares the dataset to have C. difficile colonization levels on day 5 as the values that will be predicted by regression. The features are Otu abundances at Day -1, 0, and 1 

######################## DATA PREPARATION #############################
# Features: 16S rRNA gene sequences(OTUs) in the stool
# Labels: - C. difficile CFU levels in the stool 5 days post-challenge

# Read in metadata and select only sample id, day, and cfu_d5 columns
source("code/functions.R")
metadata <- metadata %>% select(id, day, cfu_d5)
# Read in OTU table and remove label and numOtus columns
shared <- read.delim('data/process/vendors.subsample.shared', header=T, sep='\t') %>%
  select(-label, -numOtus)
# Merge metadata and OTU table.
# Filter merged dataframe to just timepoints we want to use community structure to predict C. difficile CFU 5 days post-challenge
# Then remove the sample ID and day columns
select_timepoint <- function(timepoint){
  data <- inner_join(metadata, shared, by=c("id"="Group")) %>% 
    filter(day == timepoint) %>% 
    select(-id, -day) %>%
    drop_na() %>%
    select(cfu_d5, everything()) %>%
    rename(regress=cfu_d5) %>% #rename column, the cfu_d5 values are what we want to predict
    select(-Otu0020) %>% #remove C. difficile (Otu0020) from the input data since C. difficile colonization status at day 7 is what we're predicting
    write_csv(paste0("data/process/regress_input_day", timepoint, "_data.csv"))
}
data_dn1 <- select_timepoint(-1)
data_d0 <- select_timepoint(0)
data_d1 <- select_timepoint(1)

###################################################################
