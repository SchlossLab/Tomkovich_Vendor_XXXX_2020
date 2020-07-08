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
library(ggforce) #Use geom_circle to make venn diagrams
library(gganimate) #Use to create PCoA animation
library(writexl) #For writing supplemental tables
library(ggtext)
library(glue)

### Load in metadata and make sure experiment and vendor columns are treated as factors
metadata <- read_csv("data/process/vendor_metadata.csv") %>%
  mutate(experiment=factor(experiment, levels=c("1", "2")), # Make sure experiment is treated as a factor
         cage=factor(cage, levels=c("1", "2")), #Make sure cage is treated as a factor
         vendor=factor(vendor, levels=c("Schloss", "Young", "Jackson", "Charles River", "Taconic", "Envigo"))) %>%
  filter(!grepl("Mock", id), #Remove all 5 Mock controls
         !grepl("Water", id)) %>% #Remove all 5 water controls
  #Create a variable for mouse cage (has to combine exp., vendor, and cage)
  unite(col = unique_cage, c(experiment, vendor, cage), remove = FALSE)

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

# Function to create columns for C. difficile CFU at Days 1-9
cfu_dx <- function(timepoint){
  metadata %>%
    filter(day == timepoint) %>%
    mutate("col_name" = cfu) %>%
    select(mouse_id, "col_name")
}
cfu_d1 <- cfu_dx(1) %>%
  rename(cfu_d1 = col_name)
cfu_d2 <- cfu_dx(2) %>%
  rename(cfu_d2 = col_name)
cfu_d3 <- cfu_dx(3) %>%
  rename(cfu_d3 = col_name)
cfu_d4 <- cfu_dx(4) %>%
  rename(cfu_d4 = col_name)
cfu_d5 <- cfu_dx(5) %>%
  rename(cfu_d5 = col_name) #43 values with 6 NAs
cfu_d6 <- cfu_dx(6) %>%
  rename(cfu_d6 = col_name) #34 values with 15 NAs
cfu_d7 <- cfu_dx(7) %>%
  rename(cfu_d7 = col_name) #40 values with 9 NAs
cfu_d8 <- cfu_dx(8) %>%
  rename(cfu_d8 = col_name)
cfu_d9 <- cfu_dx(9) %>%
  rename(cfu_d9 = col_name)

# add cfu_d1-9 columns to metadata----
metadata <- full_join(metadata, cfu_d1, by = "mouse_id")
metadata <- full_join(metadata, cfu_d2, by = "mouse_id")
metadata <- full_join(metadata, cfu_d3, by = "mouse_id")
metadata <- full_join(metadata, cfu_d4, by = "mouse_id")
metadata <- full_join(metadata, cfu_d5, by = "mouse_id")
metadata <- full_join(metadata, cfu_d6, by = "mouse_id")
metadata <- full_join(metadata, cfu_d7, by = "mouse_id") %>%
  #Add a column denoting C. difficile clearance status at Day 7
  mutate(clearance_status_d7 = case_when(cfu_d7 > 0 ~ "colonized",
                                         cfu_d7 == 0 ~ "not_detectable",
                                         cfu_d7 == NA ~ "NA"))
metadata <- full_join(metadata, cfu_d8, by = "mouse_id")
metadata <- full_join(metadata, cfu_d9, by = "mouse_id")

#Create a column to denote baseline_weight and weight_change from baseline_weight for each mouse----
#Calculated based on the weight recorded on D-1 of the experiment
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
color_scheme <- c("#1f78b4", "#e6ab02", "#323232", "#e7298a", "#7570b3", "#1b9e77") #Adapted from http://colorbrewer2.org/#type=qualitative&scheme=Dark2&n=6. Previous color for Jackson mice: #d95f02
color_vendors <- levels(metadata$vendor)

#Define shape scheme to differentiate experiment 1 and 2----
shape_scheme <- c(19, 17)
shape_experiment <- c(1, 2)

#Export processed metadata as a .tsv to read in to manuscript.Rmd for pulling out
#sample sizes to add to figure legends.
write_tsv(metadata, "data/process/processed_metadata.tsv")

### Calculate Standard Error----
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

#Function to find which significant otus/genera/families are shared
intersect_all <- function(a,b,...){
  Reduce(intersect, list(a,b,...))
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

#Function to calculate the median cfu values from a dataframe (x)
get_cfu_median_vendor <- function(x){
  x %>%
    group_by(vendor) %>%
    summarize(median=median(cfu)) %>%
    spread(key=vendor, value=median)
}

#Function to calculate the median weight_change values from a dataframe (x)
get_weight_median_vendor <- function(x){
  x %>%
    group_by(vendor) %>%
    summarize(median=median(weight_change)) %>%
    spread(key=vendor, value=median)
}

#Function to calculate the median agg_rel_abund values from a dataframe (x) grouped by vendor
get_rel_abund_median_vendor <- function(x){
  x %>%
    group_by(vendor) %>%
    summarize(median=median(agg_rel_abund)) %>%
    spread(key=vendor, value=median)
}

#Function to calculate the median agg_rel_abund values from a dataframe (x) grouped by day
get_rel_abund_median_day <- function(x){
  x %>%
    group_by(day) %>%
    summarize(median=median(agg_rel_abund)) %>%
    spread(key=day, value=median)
}

#Function to calculate the median agg_rel_abund values from a dataframe (x) grouped by day 7 C. diff colonization status
get_rel_abund_median_d7status <- function(x){
  x %>%
    group_by(clearance_status_d7) %>%
    summarize(median=median(agg_rel_abund)) %>%
    spread(key=clearance_status_d7, value=median)
}

#Function to calculate the median agg_rel_abund values from a dataframe (x) grouped by experiment
get_rel_abund_median_experiment <- function(x){
  x %>%
    group_by(experiment) %>%
    summarize(median=median(agg_rel_abund)) %>%
    spread(key=experiment, value=median)
}
#Function to tidy pairwise comparisons to use for adding stats to plots----
tidy_pairwise <- function(spread_pairwise){
  spread_pairwise %>%
    pivot_longer(-day, names_to = "compare", values_to = "p.adj") %>%
    separate(col = compare, c("group1", "group2"), sep = "-", remove = TRUE)
}

#Function to tidy pairwise comparisons to use for adding stats to otu plots----
tidy_pairwise_otu <- function(spread_pairwise){
  spread_pairwise %>%
    pivot_longer(-otu, names_to = "compare", values_to = "p.adj") %>%
    separate(col = compare, c("group1", "group2"), sep = "-", remove = TRUE)
}

#Function to calculate the median shannon values from a dataframe (x) grouped by vendor
get_shannon_median_vendor <- function(x){
  x %>%
    group_by(vendor) %>%
    summarize(median=median(shannon)) %>%
    spread(key=vendor, value=median)
}

#Function to calculate the median sobs (richness) values from a dataframe (x) grouped by vendor
get_sobs_median_vendor <- function(x){
  x %>%
    group_by(vendor) %>%
    summarize(median=median(sobs)) %>%
    spread(key=vendor, value=median)
}
