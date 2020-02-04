library(tidyverse)

#Load in vendors.files.csv and fix group ids by removing dashes. 
#This ensures there won't be problems running the sequencing analysis with mothur without having to change the raw data files.

#First, open up vendors.files in excel and save as .csv file

#Read in vendors.files.csv
vendor.files <- read_csv("data/raw/vendors.files.csv", col_names = FALSE) #No column names in the .files format

#Remove dashes from X1 column (the 1st column)
vendor.files <- vendor.files %>% 
  mutate(fix_id = gsub("-", "", vendor.files$X1))

#Correct id typos that were made in the MiSeq plate map, so that they will match to the metadata ids
vendor.files <- vendor.files %>% 
  mutate(fix_id = replace(fix_id, fix_id == "Y22E6E2", "Y22D6E2"), #Typo, E6 should be D6 (D stands for Day)
         fix_id = replace(fix_id, fix_id == "E212Dn1E2","E12Dn1E2")) #Josh checked his notes from loading the DNA extraction bead plates and confirmed.

#Merge any duplicate sequences by giving them the same id (sample added to 2 different wells and was sequenced twice within the same run)
vendor.files <- vendor.files %>% 
  mutate(fix_id = replace(fix_id, fix_id == "C21D0E2no2", "C21D0E2"),
         fix_id = replace(fix_id, fix_id == "E21D1E2No1", "E21D1E2"),
         fix_id = replace(fix_id, fix_id == "E21D1E2No2", "E21D1E2"),
         fix_id = replace(fix_id, fix_id == "E21Dn1E1no2", "E21Dn1E1"),
         fix_id = replace(fix_id, fix_id == "S22D1E2No1", "S22D1E2"),
         fix_id = replace(fix_id, fix_id == "S22D1E2No2", "S22D1E2"),
         fix_id = replace(fix_id, fix_id == "T12D8E1No1", "T12D8E1"),
         fix_id = replace(fix_id, fix_id == "T12D8E1No2", "T12D8E1"))

#Export .csv file
write.csv(vendor.files, file= "data/raw/vendors.files.fixed.csv")

#Open up vendors.files.fixed.csv and vendors.files
#Paste fix_id column (select all 439 samples, starting at row 2) over the 1st column with dashes in the grop name in vendors.files
#Save vendors.files now that the sample ids have been corrected by removing the dashes

#Continue running the data through mothur by moving on to code/get_good_seq.batch
