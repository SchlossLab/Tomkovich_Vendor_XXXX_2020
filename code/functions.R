library(tidyverse)
library(broom)
library(cowplot)
library(magick)
library(vegan)

### Load in metadata & create new column to give unique mouse_id based on exp. #, vendor, cage #, and mouse #.
metadata <- read_csv("data/process/vendor_metadata.csv") %>% 
  mutate(experiment=factor(experiment, levels=c("1", "2")), # Make sure experiment is treated as a factor
         vendor=factor(vendor, levels=c("Schloss", "Young", "Jackson", "Charles River", "Taconic", "Envigo"))) %>% 
  filter(!grepl("Mock", id), #Remove all 5 Mock controls
         !grepl("Water", id)) #Remove all 5 water controls

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

