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
# C21D0E2no2 (duplicate of C21D0E2?, which is present in shared_sample_names (both part of run_1, plate_3) or C21D0E1 which isn't present)
# E212Dn1E2 (duplicate of E21Dn1E2?, which is present in shared_sample_names or E12Dn1E2, which is absent in shared_sample_names? E22Dn1E2, E12Dn1E1, E11Dn1E2 is also present in shared_sample names).
# E21D1E2No1 and E21D1E2No2. Sample sequenced twice? Both were part of the same plate on run_1, plate_2.
# E21Dn1E1no2 (duplicate of E21Dn1E1?, which is present in shared_sample_names). Both a part of run_1, plate_3 and E21Dn1E2 is also present on run_1, plate_3
# S22D1E2No1 and S22D1E2No2. Sample sequenced twice? Which one to pick...Both a part of run_2, plate_1 and S22D1E2 was also sequenced on run_2, plate_3
# T12D8E1No1 and T12D8E1No2. Sample sequenced twice? Which one to pick...Both  a part of run_2, plate_3, T12D8E2 is in shared_sample_names but not run_plate_info? 

#Okay to lose T12D13E1. The 1st experiment, we collected 1 last stool sample at D13 from this mouse since it was the only mouse still colonized.
# Y13D8E2. Okay to lose this sample? There was only a Y13 mouse in the 1st experiment and Y13D8E1 is already listed in shared_sample_names. 
#Looking back at the plate layout maps Y13D9E1 (Y13_D9_E1) was listed twice, so maybe this sample was really a duplicate of that one.
#Y13D9E1 shows up in run_plate_info duplicates, while Y13D8E2 is not present in run_plate_info

#Check for any duplicate values in final metadata data frame:
duplicated <- metadata %>% 
  filter(duplicated(id))
# One sample was listed twice: Y21D0E2. 2nd instance should actually be  Y21D1E2. 
# Corrected 2nd instance of Y21D0E2 to Y21D1E2 in the metadata.v2.csv

#Load in ids, MiSeq run info, plate #, and plate locations for the sequenced samples----
run_plate_info <- read_csv("data/process/vendor_16S_run_plate_info.csv") %>% 
  mutate(id = gsub("_", "", original_id)) %>%  # remove underscores from plate_map id labels to match sample ids used for mothur pipeline
  filter(!str_detect(id, "Mock|Water|NoSample")) # remove rows with Mock, Water, or NoSample ids
#455 samples based on the plate map files for the 2 MiSeq runs that sequenced vendor project samples

#Check for duplicated ids based on labels from the plate map
duplicated <- run_plate_info %>% 
  filter(duplicated(id)) %>% 
  pull(id)
# 23 duplicates
#List of duplicates:  [1] "J12D9E1" "S11D9E1" "C12D9E1" "J22D9E1" "S21D9E1" "J11D9E1"
#[7] "C21D9E1" "T11D9E1" "S22D9E1" "E22D9E1" "E21D9E1" "Y13D9E1"
#[13] "T21D9E1" "S12D9E1" "J21D9E1" "C11D9E1" "Y12D9E1" "T12D9E1"
#[19] "Y21D9E1" "E12D9E1" "C22D9E1" "E11D9E1" "Y22D9E1"

#Pull out duplicated samples in run_plate_info
run_plate_info_duplicates <-  run_plate_info %>% 
  filter(id %in% c(duplicated)) 
#Returns 46 samples
run1_duplicates <- run_plate_info_duplicates %>% filter(run == "run_1")  #23 samples
run2_duplicates <- run_plate_info_duplicates %>% filter(run == "run_2")  #23 samples
inner_join(run1_duplicates, run2_duplicates, by = "id") #23 samples
#So all 23 of these samples were sequenced twice, once on the 1st MiSeq_run, and a 2nd time on the 2nd MiSeq_Run


#Pull out samples with no in name (also represent duplicates)
run_plate_info_no1or2 <- run_plate_info %>% 
  filter(str_detect(id, "No|no"))
#Returns 8 samples
# C21D0E2no2 & E21D1E2No2 are the only pair listed without a corresponding no1

#Check what's shared or different between mothur outputted sample ids in the .shared files and ids from the plate map----
matching_platemap <- inner_join(run_plate_info, shared_sample_names, by = c('id' = 'Group'))
#433 samples match between plate_map & shared file.
matching_platemap_duplicates <- matching_platemap %>%  filter(duplicated(id)) 
#23 of these matches are duplicates

#What's missing from plate_map compared to shared_sample_names----
missing_platemap <- anti_join(shared_sample_names, run_plate_info, by = c('Group' = "id")) %>% pull(Group)
#19 samples: "C11D8E2" "C21D8E2" "C22D8E2" "E12D8E2" "E21D8E2" "E22D8E2" "J11D8E2" "J21D8E2" "J22D8E2" "S12D8E2" "S21D8E2" "S22D8E2" "T11D8E2" "T12D8E2" "T21D8E2" "Y12D8E2" "Y13D8E2" "Y21D8E2" "Y22D8E2" 
#All have corresponding .fastq files in the raw data folder but not a corresponding pair of .fastq files in the BaseCalls folder of either sequencing run

#What's missing from shared_sample_names compared to plate_map----
missing_shared_samples <- anti_join(run_plate_info, shared_sample_names, by = c('id' = 'Group')) %>% pull(id)
#20 samples: "C21Dn1E1" "S21Dn1E2" "C11D0E2" "E21D0E1" "Y11D5E2" "C11D2E2" "C12D2E2" "T11D7E2" "Y21D6E2" "T21D4E2" "T22D7E2" "C22D7E2" "C22D2E2" "Y21D1E2" "Y22D4E2" "Y22D3E2" "S22D3E2ANDT22D3E2" "S21D7E1" "T21D7E1" "S22D3E1"
#All but "S22D3E2ANDT22D3E2" are found as .fastq files in raw data folder. Suspect these samples were dropped during subsampling?

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