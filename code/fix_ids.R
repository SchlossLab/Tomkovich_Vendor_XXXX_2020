library(tidyverse)

#Load in vendors.files.csv and fix group ids by removing dashes. 
#This ensures there won't be problems running the sequencing analysis with mothur without having to change the raw data files.

#First, open up vendors.files in excel and save as .csv file

#Read in vendors.files.csv
vendor.files <- read_csv("data/raw/vendors.files.csv", col_names = FALSE) #No column names in the .files format

#Remove dashes from X1 column (the 1st column)
vendor.files <- vendor.files %>% 
  mutate(fix_id = gsub("-", "", vendor.files$X1))

#Export .csv file
write.csv(vendor.files, file= "data/raw/vendors.files.fixed.csv")

#Open up vendors.files.fixed.csv and vendors.files
#Paste fix_id column (select all 439 samples, starting at row 2) over the 1st column with dashes in the grop name in vendors.files
#Save vendors.files now that the sample ids have been corrected by removing the dashes

#Continue running the data through mothur by moving on to code/get_good_seq.batch
