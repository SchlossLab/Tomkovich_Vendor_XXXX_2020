#Splits the sequenced samples into individual timepoints or vendors. 
#Used to fill in groups parameter of dX_dist_PCoA.batch and dX_vendor_PCoA.batch scripts. 

source("code/functions.R")

#Full list of samples with sequence data that includes all timepoints
pcoa_data <- read_tsv("data/process/vendors.subsample.thetayc.ave.pcoa.axes") %>%
  select(group, axis1, axis2) %>% #Limit to 2 PCoA axes
  rename(id = group) %>% #group is the same as id in the metadata data frame
  right_join(metadata, by= "id") %>% #merge metadata and PCoA data frames
  filter(!is.na(axis1)) #Remove all samples that weren't sequenced or were sequenced and didn't make the subsampling cutoff

#Day -1 samples----
pcoa_dn1 <- pcoa_data %>% 
  filter(day == -1) %>% 
  pull(id) %>% 
  noquote() #Remove quotations from the characters

#Concatenate output and add - between each sample.
dn1_ids <- (paste(pcoa_dn1, collapse = "-"))
#Day-1 ids:
#Output with the quotes removed, which is what we'll use for groups argument for the dist.shared command in mothur
#S11Dn1E2-S12Dn1E2-S22Dn1E2-Y11Dn1E2-Y12Dn1E2-Y21Dn1E2-Y22Dn1E2-J12Dn1E2-J21Dn1E2-J22Dn1E2-C11Dn1E2-C12Dn1E2-C21Dn1E2-C22Dn1E2-T11Dn1E2-T12Dn1E2-T21Dn1E2-T22Dn1E2-E11Dn1E2-E12Dn1E2-E21Dn1E2-E22Dn1E2-S11Dn1E1-S12Dn1E1-S21Dn1E1-Y11Dn1E1-Y12Dn1E1-Y21Dn1E1-Y22Dn1E1-J11Dn1E1-J12Dn1E1-C12Dn1E1-C22Dn1E1-T11Dn1E1-T12Dn1E1-T21Dn1E1-T22Dn1E1-E12Dn1E1-E21Dn1E1-E22Dn1E1

#Dn1 Schloss samples----
pcoa_dn1_s <- pcoa_data %>% 
  filter(day == -1 & vendor == "Schloss") %>% 
  pull(id) %>% 
  noquote() 
dn1_s_ids <- (paste(pcoa_dn1_s, collapse = "-"))
#Day -1 Schloss ids:
#S11Dn1E2-S12Dn1E2-S22Dn1E2-S11Dn1E1-S12Dn1E1-S21Dn1E1

#Dn1 Young samples----
pcoa_dn1_y <- pcoa_data %>% 
  filter(day == -1 & vendor == "Young") %>% 
  pull(id) %>% 
  noquote() 
dn1_y_ids <- (paste(pcoa_dn1_y, collapse = "-"))
#Day -1 Young ids:
#Y11Dn1E2-Y12Dn1E2-Y21Dn1E2-Y22Dn1E2-Y11Dn1E1-Y12Dn1E1-Y21Dn1E1-Y22Dn1E1

#Dn1 Jackson samples----
pcoa_dn1_j <- pcoa_data %>% 
  filter(day == -1 & vendor == "Jackson") %>% 
  pull(id) %>% 
  noquote() 
dn1_j_ids <- (paste(pcoa_dn1_j, collapse = "-"))
#Day -1 Jackson ids:
#J12Dn1E2-J21Dn1E2-J22Dn1E2-J11Dn1E1-J12Dn1E1

#Dn1 Charles River samples----
pcoa_dn1_c <- pcoa_data %>% 
  filter(day == -1 & vendor == "Charles River") %>% 
  pull(id) %>% 
  noquote() 
dn1_c_ids <- (paste(pcoa_dn1_c, collapse = "-"))
#Day -1 Charles River ids:
#C11Dn1E2-C12Dn1E2-C21Dn1E2-C22Dn1E2-C12Dn1E1-C22Dn1E1

#Dn1 Taconic samples----
pcoa_dn1_t <- pcoa_data %>% 
  filter(day == -1 & vendor == "Taconic") %>% 
  pull(id) %>% 
  noquote() 
dn1_t_ids <- (paste(pcoa_dn1_t, collapse = "-"))
#Day -1 Taconic ids:
#T11Dn1E2-T12Dn1E2-T21Dn1E2-T22Dn1E2-T11Dn1E1-T12Dn1E1-T21Dn1E1-T22Dn1E1

#Dn1 Envigo samples----
pcoa_dn1_e <- pcoa_data %>% 
  filter(day == -1 & vendor == "Envigo") %>% 
  pull(id) %>% 
  noquote() 
dn1_e_ids <- (paste(pcoa_dn1_e, collapse = "-"))
#Day -1 Envigo ids:
#E11Dn1E2-E12Dn1E2-E21Dn1E2-E22Dn1E2-E12Dn1E1-E21Dn1E1-E22Dn1E1

#Day 0 samples----
pcoa_d0 <- pcoa_data %>% 
  filter(day == 0) %>% 
  pull(id) %>% 
  noquote() #Remove quotations from the characters

#Concatenate output and add - between each sample.
d0_ids <- (paste(pcoa_d0, collapse = "-"))
#Day0 ids:
#Output with the quotes removed, which is what we'll use for groups argument for the dist.shared command in mothur
#S11D0E2-S12D0E2-S21D0E2-S22D0E2-Y11D0E2-Y12D0E2-Y21D0E2-Y22D0E2-J11D0E2-J12D0E2-J21D0E2-J22D0E2-C21D0E2-C22D0E2-T11D0E2-T12D0E2-T22D0E2-E11D0E2-E12D0E2-E22D0E2-S12D0E1-S21D0E1-S22D0E1-Y12D0E1-Y13D0E1-Y21D0E1-Y22D0E1-J11D0E1-J12D0E1-J21D0E1-J22D0E1-C11D0E1-C12D0E1-C22D0E1-T11D0E1-T21D0E1-E11D0E1-E12D0E1-E22D0E1

#Day 1 samples----
pcoa_d1 <- pcoa_data %>% 
  filter(day == 1) %>% 
  pull(id) %>% 
  noquote() #Remove quotations from the characters

#Concatenate output and add - between each sample.
d1_ids <- (paste(pcoa_d1, collapse = "-"))
#Day1 ids:
#Output with the quotes removed, which is what we'll use for groups argument for the dist.shared command in mothur
#S11D1E2-S12D1E2-S21D1E2-S22D1E2-Y11D1E2-Y12D1E2-J11D1E2-J12D1E2-J22D1E2-C11D1E2-C12D1E2-C21D1E2-C22D1E2-T11D1E2-T21D1E2-T22D1E2-E11D1E2-E12D1E2-E21D1E2-E22D1E2-S11D1E1-S12D1E1-S21D1E1-S22D1E1-Y11D1E1-Y12D1E1-Y13D1E1-Y21D1E1-Y22D1E1-J21D1E1-J22D1E1-C11D1E1-C12D1E1-C21D1E1-C22D1E1-T11D1E1-T12D1E1-T21D1E1-E11D1E1-E12D1E1-E21D1E1-E22D1E1
