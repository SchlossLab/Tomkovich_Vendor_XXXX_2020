library(tidyverse)
library(broom)
library(cowplot)
library(magick)
library(vegan)
library(reshape2)
library(knitr)
library(rmarkdown) 
library(gtools)
library(ggpubr)

### Load in metadata & create new column to give unique mouse_id based on exp. #, vendor, cage #, and mouse #.
metadata <- read_csv("data/process/vendor_metadata.csv") %>% 
  mutate(experiment=factor(experiment, levels=c("1", "2")), # Make sure experiment is treated as a factor
         vendor=factor(vendor, levels=c("Schloss", "Young", "Jackson", "Charles River", "Taconic", "Envigo"))) %>% 
  filter(!grepl("Mock", id), #Remove all 5 Mock controls
         !grepl("Water", id)) #Remove all 5 water controls

#Quantify C. diff cfu based on the colonies counted from the 2 different dilutions plated. Formula based on protocol used by the Young lab----
#0s should only be kept when the -1 dilution was checked, as that represents the limit of detection.
# If -1 dilution was not plated, we can not say for sure whether C. diff CFU for a mouse is really zero. Thus, zeros for any dilutions greater than the -1 dilution were converted to NAs
metadata <- metadata %>% 
  mutate(cfu1 = count1 * 20 * 1 / (10 ^ dilution1), cfu2 = count2 * 20 * (1 / (10 ^ dilution2))) %>% # Quantify CFU/g for each dilution that was plated
  select(-starts_with("count")) # gets rid of count columns since these are now represented as cfu1 and cfu2 columns

#Make sure 0s are true 0s, meaning the -1 dilution (which is the limit of detection) was plated. Any 0s for dilutions above the -1 dilution should be transformed to NAs. 
#Number of 0s in cfu1 should equal number of -1s in dilution1. Since -1 dilution was never plated in dilution2 column, don't have to worry about counting any zeros from that column.
#Quantify how many instances we have 0s for cfu1, cfu2, etc. so we can be sure we transformed data correctly in the next step.
cfu1_0s <- metadata %>% filter(cfu1 == 0) #229/563 instances
cfu1_d <- metadata %>% filter(dilution1 == -1 & cfu1 == 0) #184 instances, meaning there should only be 184 0s.
cfu2_0s <- metadata %>% filter(cfu2 == 0) #133/563 instances. Most of these should be transformed to NAs
cfu_nas <- map(metadata, ~sum(is.na(.))) #151 for cfu1, 234 for cfu2

#Transform data so that 0s from -1 dilution remain 0s, but 0s for any dilutions beyond -1 become NA.----
metadata <- metadata %>%
  mutate(cfu1 = ifelse(cfu1 == 0 & dilution1 != -1, NA, cfu1)) %>% #Keeps 0s for -1 dilution, replaces 0s from any other dilution with NA
  mutate_at(vars(cfu2), ~replace(., . == 0, NA)) %>% #Changes all 0s for cfu2 to NA since the -1 dilution (the limit of detection) was not plated.
  group_by(id, experiment, mouse_id, day) %>% 
  mutate(cfu = mean(c(cfu1, cfu2), na.rm = TRUE)) %>% #Create a final cfu per ID (combination of mouse ID & date sample was collected) based on the average CFU/g for cfu1 and cfu2
  mutate(cfu = na_if(cfu, "NaN")) %>% #Changes the NaNs in cfu column back to Nas
  ungroup()
cfu1_0s <- metadata %>% filter(cfu1 == 0) #Now only 184 instances, which is what we predicted on line 21
cfu2_0s <- metadata %>% filter(cfu2 == 0) #0 instances of 0.
cfu_nas_final <- map(metadata, ~sum(is.na(.))) #196 for cfu1, 367 instances for cfu2. #161 for cfu

# Create columns for C. difficile CFU at Day 5 and Day 6 since these were the 2 timepoints where there were significant differences across vendors and the most significant pairwise.wilcox values (6 pairs)
cfu_d5 <- metadata %>% 
  filter(day == 5) %>% 
  mutate(cfu_d5 = cfu) %>% 
  select(mouse_id, cfu_d5) #43 values with 6 NAs
cfu_d6 <- metadata %>% 
  filter(day == 6) %>% 
  mutate(cfu_d6 = cfu) %>% 
  select(mouse_id, cfu_d6) #34 values with 15 NAs
cfu_d7 <- metadata %>% 
  filter(day == 7) %>% 
  mutate(cfu_d7 = cfu) %>% 
  select(mouse_id, cfu_d7) #40 values with 5 NAs. 17 mice with 0, 23 still colonized. Close to half have cleared. By D8 33 mice with 0s, 3 still colonized, 13 NAs

# add cfu_d5, d6, d7 columns to metadata----
metadata <- full_join(metadata, cfu_d5, by = "mouse_id")
metadata <- full_join(metadata, cfu_d6, by = "mouse_id")
metadata <- full_join(metadata, cfu_d7, by = "mouse_id") %>% 
  #Add a column denoting C. difficile clearance status at Day 7
  mutate(clearance_status_d7 = case_when(cfu_d7 > 0 ~ "colonized",
                                         cfu_d7 == 0 ~ "cleared",
                                         cfu_d7 == NA ~ "NA"))

#Create a column to denote baseline_weight and weight_change from baseline_weight for each mouse----
#Percent baseline weight for each mouse will be calculated 
#based on the weight recorded on D-1 of the experiment
baseline_weight <- metadata %>% select(mouse_id, weight, day) %>% 
  filter(day == -1) %>% 
  mutate(baseline_weight = weight) %>% 
  select(mouse_id, baseline_weight)
#Join baseline data frame to main metadata
metadata <- inner_join(metadata, baseline_weight, by = "mouse_id") %>% 
  #Calculate weight change(g) for each mouse based on D-1 weight.
  group_by(mouse_id, day) %>% 
  mutate(weight_change = weight-baseline_weight) %>% 
  ungroup()
  
#Check to make sure sample ids on shared file match the sample ids in the metadata file----
#Figure out number of samples that were sequenced:
shared_sample_names <- read.table('data/process/vendors.subsample.shared', 
                                  sep = '\t', header = T, stringsAsFactors = F) %>% 
  select(Group) # select column with sample names of shared file which we will need to match our metadata to
#404 samples with sequence data. 

#Figure out which samples match/don't match between the shared and metadata files
matching <- inner_join(metadata, shared_sample_names, by = c('id' = 'Group'))
#404 samples match between metadata & shared file.

#Check for any duplicate ids in final metadata data frame:
duplicated <- metadata %>% 
  filter(duplicated(id))
#0 duplicated samples in metadata

#Check for any duplicate ids in shared_sample_names data frame:
duplicated <- duplicated(shared_sample_names)
#0 duplicated samples in shared_sample_names

#Dataframe of just the samples that were sequenced on the MiSeq
sequenced <- metadata %>% 
  filter(!is.na(run)) 
#423 samples were sequenced on the MiSeq

dropped_by_miseq <- nrow(sequenced)-nrow(matching)
#19 samples were dropped during subsampling to 5437 sequences

#IDs of the 19 samples that were dropped after subsampling to 5437 sequences
dropped_by_miseq_ids <- anti_join(sequenced, matching, by = "id") %>% pull(id)


#Define color scheme----
color_scheme <- c("#1f78b4", "#e6ab02", "#d95f02", "#e7298a", "#7570b3", "#1b9e77") #Adapted from http://colorbrewer2.org/#type=qualitative&scheme=Dark2&n=6
color_vendors <- levels(metadata$vendor)

### Caluclate Standard Error----
calc_se <- function(stdev, n){
  stdev/sqrt(n)
}

### Calculate 95% confidence intervals----
lower_ci <- function(mean, se, n, conf_level = 0.95){
  lower_ci <- mean - qt(1 - ((1 - conf_level) / 2), n - 1) * se
}

upper_ci <- function(mean, se, n, conf_level = 0.95){
  upper_ci <- mean + qt(1 - ((1 - conf_level) / 2), n - 1) * se
}

#Function to have y-axis in scientific notation----
fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "%*%10^", l)
  # return this as an expression
  parse(text=l)
}

#Functions for machine learning analysis from Begum's paper:
#Source: https://github.com/SchlossLab/Topcuoglu_ML_XXX_2019/blob/master/code/learning/functions.R
# -------------------- Make performance files tidy------------------>
# Instead of 2 columns with names cv_aucs and test_aucs
# We will have 1 column with name Performance that tells us if test or cv
melt_data <-  function(data) {
  data_melt <- data %>%
    melt(measure.vars=c('cv_aucs', 'test_aucs')) %>%
    rename(AUC=value) %>%
    mutate(Performance = case_when(variable == "cv_aucs" ~ 'Cross-validation', variable == "test_aucs" ~ 'Testing')) %>%
    group_by(Performance)
  return(data_melt)
}
# -------------------------------------------------------------------->


# -------------------- Read files ------------------------------------>
# Read in files as delim that are saved in a list with a pattern
read_files <- function(filenames){
  for(file in filenames){
    # Read the files generated in main.R
    # These files have cvAUCs and testAUCs for 100 data-splits
    data <- read.delim(file, header=T, sep=',')
  }
  return(data)
}
# -------------------------------------------------------------------->

#Function to calculate the mean cfu values from a dataframe (x) 
get_cfu_mean_vendor <- function(x){
  x %>%
    group_by(vendor) %>%
    summarize(mean=mean(cfu)) %>%
    spread(key=vendor, value=mean)
}

