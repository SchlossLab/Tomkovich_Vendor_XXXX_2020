#Source: https://github.com/SchlossLab/ML_pipeline_microbiome/blob/master/test/code/prepare_input_data.R

#Script prepares the dataset to have C. difficile colonization status on day 7 post-challenge as the outcome. The features are family abundances at Day -1, 0, and 1.
#The C. difficile Otu (Otu 20) was filtered out prior to creating input data for relative abundances at the family level

######################## DATA PREPARATION #############################
# Features: 16S rRNA gene sequences at the family level in the stool on Day -1 of the experiment
# Labels: - C. difficile colonization status in the stool 7 days post-challenge

# Read in metadata and select only sample id, day, and clearance_status_d7 columns
source("code/functions.R")
metadata <- metadata %>% select(id, day, clearance_status_d7)
# Read in OTU table and remove label and numOtus columns
shared <- read.delim('data/process/vendors.subsample.shared', header=T, sep='\t') %>%
  select(-label, -numOtus) %>% 
  select(-Otu0020) %>% #remove C. difficile (Otu0020) from the input data since C. difficile colonization status at day 7 is what we're predicting
  gather(-Group, key=OTU, value=count)

taxonomy <- read_tsv(file="data/process/vendors.taxonomy") %>%
  select(-Size) %>%
  mutate(Taxonomy=str_replace_all(Taxonomy, "\\(\\d*\\)", "")) %>%
  mutate(Taxonomy=str_replace_all(Taxonomy, ";$", "")) %>%
  separate(Taxonomy, c("kingdom", "phylum", "class", "order", "family", "genus"), sep=';') %>%
  select(OTU, "family") %>%
  rename(taxon = family)

unique_taxonomy <- taxonomy %>%
  select(taxon) %>%
  unique() %>%
  mutate(otu = paste0("Otu", str_pad(1:nrow(.), width=nchar(nrow(.)), pad="0")))

family_shared <- inner_join(shared, taxonomy, by="OTU") %>%
  group_by(taxon, Group) %>%
  summarize(count = sum(count)) %>%
  ungroup() %>%
  inner_join(., unique_taxonomy) %>%
  select(-taxon) %>%
  spread(otu, count) %>%
  mutate(label="family", numOtus=ncol(.)-1) %>%
  select(label, Group, numOtus, everything())
write_tsv(family_shared, path = "data/process/vendors.subsample.family.shared")

select(family_shared, -label, -numOtus) %>%
  gather(otu, count, -Group) %>%
  group_by(otu) %>%
  summarize(count=sum(count)) %>%
  inner_join(., unique_taxonomy) %>%
  rename("OTU"="otu", "Size"="count", "Taxonomy"="taxon") %>%
  write_tsv(path ="data/process/vendors.family.taxonomy")

#Modify family_shared file so that only Group and Otu columns remain:
family_shared <- family_shared %>% 
  select(-label, -numOtus)

# Function that will do the following
# Merge metadata and shared file with taxa classified at the family level.
# Filter merged dataframe to just timepoints we want to use community structure to predict C. difficile CFU 7 days post-challenge
# Then remove the sample ID and day columns
select_timepoint <- function(timepoint){
  data <- inner_join(metadata, family_shared, by=c("id"="Group")) %>% 
    filter(day == timepoint) %>% 
    select(-id, -day) %>%
    drop_na() %>%
    select(clearance_status_d7, everything()) %>%
    rename(dx=clearance_status_d7) %>% #rename column, the clearance_status_d7 values are what we want to predict

    write_csv(paste0("data/process/classification_input_day", timepoint, "_data_family.csv"))
}
class_data_dn1 <- select_timepoint(-1)
class_data_d1 <- select_timepoint(1)
class_data_d0 <- select_timepoint(0)

###################################################################
