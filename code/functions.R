library(tidyverse)
library(broom)
library(cowplot)
library(magick)

### Load in metadata & create new column to give unique mouse_id based on exp. #, vendor, cage #, and mouse #.
metadata <- read_csv("data/process/vendor_metadata.v2.csv") %>% 
  unite(col = mouse_id, c(experiment, vendor, cage, mouse), remove = FALSE) %>% 
  mutate(experiment=factor(experiment, levels=c("1", "2")), # Make sure experiment is treated as a factor
         vendor=factor(vendor, levels=c("Schloss", "Young", "Jackson", "Charles River", "Taconic", "Envigo"))) 

#Check to make sure sample ids on shared file match the sample ids in the metadata file----
#Figure out number of samples that were sequenced:
shared_sample_names <- read.table('data/process/vendors.subsample.shared', 
                                  sep = '\t', header = T, stringsAsFactors = F) %>% 
  select(Group) # select column with sample names of shared file which we will need to match our metadata to
#429 samples were sequenced. 

#Figure out which samples match/don't match between the shared and metadata files
matching <- inner_join(metadata, shared_sample_names, by = c('id' = 'Group'))
#418 samples match between metadata & shared file.

# Y22E6E2 should be Y22D6E2. Correct this name on metadata in R so that it will match
# all mothur files that will be joined to metadata
metadata <- metadata %>% 
  mutate(id = replace(id, id == "Y22D6E2", "Y22E6E2"))

#What's missing from metadata compared to shared_sample_names
missing <- anti_join(shared_sample_names, metadata, by = c('Group' = "id"))
#11 samples on shared_sample_names that don't correspond to ids listed in the metadata.
#Previously 14 samples were missing, but since 2 samples were mislabeled in metadata.csv. I Updated the metadata file to fix mislabeled samples in id column.
#Checked again and number of missing samples dropped to 12 as expected. 
# Sample ids that were corrected in metadata.csv:
# 2nd istance of C11D0E2 to C11D1E2
# Update S22Dn2E2 to S22Dn1E2 (none of the samples were sequenced at Dn2, likely a typo)

#Some samples were sequenced twice. Checked with Lucas and the sample should be the same, just seqeunced across 2 wells.
# C21D0E2no2 (duplicate of C21D0E2?, which is present in shared_sample_names)
# E212Dn1E2 (duplicate of E21Dn1E2?, which is present in shared_sample_names. E22Dn1E2 is also present in shared_sample names).
# E21D1E2No1 and E21D1E2No2. Sample sequenced twice? Which one to pick...
# E21Dn1E1no2 (duplicate of E21Dn1E1?, which is present in shared_sample_names)
# S22D1E2No1 and S22D1E2No2. Sample sequenced twice? Which one to pick...
# T12D8E1No1 and T12D8E1No2. Sample sequenced twice? Which one to pick...

#Okay to lose T12D13E1. The 1st experiment, we collected 1 last stool sample at D13 from this mouse since it was the only mouse still colonized.
# Y13D8E2. Okay to lose this sample? There was only a Y13 mouse in the 1st experiment and Y13D8E1 is already listed in shared_sample_names. 
#Looking back at the plate layout maps Y13D9E1 (Y13_D9_E1) was listed twice, so maybe this sample was really a duplicate of that one.
# 13 sets of duplicate samples on plate layout sheet. Only 2 likely apply to some of the missing samples (C21D0E2 & E212Dn1E2)

#Check for any duplicate values in final metadata data frame:
duplicated <- metadata %>% 
  filter(duplicated(id))
# One sample is listed twice: Y21D0E2. 2nd instance should actually be  Y21D1E2. 
# Corrected 2nd instance of Y21D0E2 to Y21D1E2 in the metadata.csv

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