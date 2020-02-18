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
#404 samples with sequence data. 

#Figure out which samples match/don't match between the shared and metadata files
matching <- inner_join(metadata, shared_sample_names, by = c('id' = 'Group'))
#404 samples match between metadata & shared file.

# 3 mislabeled sample ids that were corrected in metadata.v2.csv:
# 2nd instance of C11D0E2 to C11D1E2
# Update S22Dn2E2 to S22Dn1E2 (none of the samples were sequenced at Dn2, likely a typo)
# 2nd instance of Y21D0E2 should actually be  Y21D1E2. 

#Check for any duplicate values in final metadata data frame:
duplicated <- metadata %>% 
  filter(duplicated(id))
#0 duplicated samples in metadata

duplicated <- duplicated(shared_sample_names)
#0 duplicated samples in shared_sample_names

#Load in ids, MiSeq run info, plate #, and plate locations for the sequenced samples----
run_plate_info <- read_csv("data/process/vendor_16S_run_plate_info.csv") %>% 
  mutate(id = gsub("_", "", original_id)) %>%  # remove underscores from plate_map id labels to match sample ids used for mothur pipeline
  filter(!str_detect(id, "Mock|Water|NoSample")) # remove rows with Mock, Water, or NoSample ids
#453 samples based on the plate map files for the 2 MiSeq runs that sequenced vendor project samples

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
#Discussed these with Josh, he thought these were D8E2 samples instead but since the mouse ids aren't a perfect match
#between experiements, we can't be confident of the label. Josh said there was not a way to differentiate exp. 1 versus exp. 2
#He also said after samples were aliquoted into plates the tubes were thrown out.
#Since we can't be confident of the identies of these sequences, I excluded them from the analysis by removing them from the data/raw/ directory

#Pull out samples with no in name (also represent duplicates)
run_plate_info_no1or2 <- run_plate_info %>% 
  filter(str_detect(id, "No|no"))
#Returns 8 samples
# C21D0E2no2 & E21D1E2No2 are the only pair listed without a corresponding no1
#All of these were corrected so that both sequences would be attributed to the same sample
#This was done with code/make_file_fix_ids.batch and fix_ids.R. vendors.files was updated as directed before mothur analysis was continued with code/get_good_seqs.batch

#Check what's shared or different between mothur outputted sample ids in the .shared files and ids from the plate map----
matching_platemap <- inner_join(run_plate_info, shared_sample_names, by = c('id' = 'Group'))
#422 samples match between plate_map & shared file.
matching_platemap_duplicates <- matching_platemap %>%  filter(duplicated(id)) 
#23 of these matches are duplicates

#What's missing from run_plate_info compared to shared_sample_names----
missing_platemap <- anti_join(shared_sample_names, run_plate_info, by = c('Group' = "id")) %>% pull(Group)
#5 samples: "E12Dn1E2" "E21D1E2" "S22D1E2" "T12D8E1" "Y22D6E2" 
#All of these were duplicates/typos that were corrected so that both sequences would be attributed to the same sample in the case of duplicates, or for typos, the sample would match up with the corresponding metadata entry
#This was done with code/make_file_fix_ids.batch and fix_ids.R. vendors.files was updated as directed before mothur analysis was continued with code/get_good_seqs.batch

#What's missing from shared_sample_names compared to run_plate_info----
missing_shared_samples <- anti_join(run_plate_info, shared_sample_names, by = c('id' = 'Group')) %>% pull(id)
#31 samples: "C21Dn1E1" "E212Dn1E2" "S21Dn1E2""E21Dn1E1no2""C21D0E2no2""C11D0E2""E21D0E1""Y11D5E2""S22D1E2No1""C11D2E2""C12D2E2""T11D7E2""Y21D6E2""T21D4E2""T22D7E2""C22D7E2""C22D2E2""S22D1E2No2""Y22E6E2""E21D1E2No2""Y21D1E2""E21D1E2No1""Y22D4E2""Y22D3E2" "S22D3E2ANDT22D3E2" "T12D13E1""S21D7E1""T21D7E1""S22D3E1""T12D8E1No2""T12D8E1No1" 
# "S22D3E2ANDT22D3E2" sequences were removed from data/raw since it represents sequences from 2 different samples.
# "T12D13E1"sequences were removed from data/raw since it represents sequences collected 13 days post-infection and was the only mouse with a sample collected at that timepoint.
# 10 Samples represent corrections made with code/make_file_fix_ids.batch and fix_ids.R. vendors.files was updated as directed before mothur analysis was continued with code/get_good_seqs.batch: "E212Dn1E2" "E21Dn1E1no2" "C21D0E2no2" "S22D1E2No1" "S22D1E2No2" "Y22E6E2" "E21D1E2No2" "E21D1E2No1" "T12D8E1No2""T12D8E1No1"
# 19 Samples dropped after rarefying to 5437 sequences: "C21Dn1E1" "S21Dn1E2" "C11D0E2" "E21D0E1" "Y11D5E2" "C11D2E2" "C12D2E2" "T11D7E2" "Y21D6E2" "T21D4E2" "T22D7E2" "C22D7E2" "C22D2E2" "Y21D1E2" "Y22D4E2" "Y22D3E2" "S21D7E1" "T21D7E1" "S22D3E1"

#Check if anything is missing from metadata compared to list of samples on run_plate_info
missing_metadata_v_runplate <- anti_join(run_plate_info, metadata, by = "id") %>% pull(id)
#12 samples that were previous identified above:  "E212Dn1E2" "E21Dn1E1no2" "C21D0E2no2" "S22D1E2No1" "S22D1E2No2" "Y22E6E2" "E21D1E2No2"       
# "E21D1E2No1" "S22D3E2ANDT22D3E2" "T12D13E1" "T12D8E1No2" "T12D8E1No1" 
#Check for the opposite situation
missing_runplate_v_metadata <- anti_join(metadata, run_plate_info, by = "id") %>% pull(id)
#145 samples are in metadata but not in the plate map. 

#Join run/plate/plate_location info to metadata file----
run_plate_info <- run_plate_info %>% 
  select(id, run, plate, plate_location)
metadata <- left_join(metadata, run_plate_info, by = "id") #Increased by 23 rows thanks to the 23 duplicates across 2 runs
duplicated_v2 <- metadata %>% 
  filter(duplicated(id)) #Returns the 23 duplicates (all from D9 timepoint)
metadata <- distinct(metadata, id, .keep_all = TRUE) %>% #Keeps just 1 of the duplicates
  write_xlsx(path = "data/process/metadata.v3.xlsx", format_headers = FALSE)


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

