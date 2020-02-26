source("code/functions.R")

all_data <- read_tsv("data/process/vendors.subsample.shared") %>%
  left_join(metadata, by = c("Group" = "id"))

all <- all_data %>%
  select(Group, starts_with("Otu")) %>%
  column_to_rownames(var = "Group")

exp1 <- all_data %>%
  filter(experiment == 1) %>%
  select(Group, starts_with("Otu")) %>%
  column_to_rownames(var = "Group")

exp2 <- all_data %>%
  filter(experiment == 2) %>%
  select(Group, starts_with("Otu")) %>%
  column_to_rownames(var = "Group")

#Function to calculate PCoA scores based on a specified dataframe of samples
#Must supply a distance matrix, otherwise vegdist will be used to find the dissimilarities
#Use avgdist() to create distance matrix using subsampling and iterations to mimic dist.shared calculator in mothur
calculate_pcoa <- function(subset_data){
  dist_matrix <- avgdist(subset_data, sample = 5347, distfun = vegdist, meanfun = mean, transf = NULL, iterations = 1000, dmethod = "bray")
  pcoa <- capscale(dist_matrix~1, distance = "bray")
  pcoa_data <- as.data.frame(scores(pcoa, display="sites")) %>% #Note scores(exp1_pcoa, display="species") #Displays scores based on Otus
    mutate("id" = rownames(.)) %>% #Convert samples as rownames to entries in id coloumn
    right_join(metadata, by= "id") %>% #merge metadata and PCoA data frames
    filter(!is.na(MDS1)) #Remove all samples that weren't sequenced or were sequenced and didn't make the subsampling cutoff
}

all_pcoa <- calculate_pcoa(all)
exp1_pcoa <- calculate_pcoa(exp1)
exp2_pcoa <- calculate_pcoa(exp2)

exp_days <- c(-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
for(d in exp_days){
  name <- paste("day", d, "_pcoa", sep = "")
  subset_data <- all_data %>% 
    filter(day == d) %>% 
    select(Group, starts_with("Otu")) %>%
    column_to_rownames(var = "Group")
  assign(name, calculate_pcoa(subset_data))
}

#Function to plot PCoA data for all vendors----
plot_pcoa <- function(df){
  ggplot(df, aes(x=MDS1, y=MDS2, color = vendor, alpha = day)) +
    geom_point(size=2) +
    scale_colour_manual(name=NULL,
                        values=color_scheme,
                        breaks=color_vendors,
                        labels=color_vendors)+
#    scale_alpha_continuous(range = c(.3, 1),
#                           breaks= c(2, 4, 6, 8, 10),
#                           labels=c(2, 4, 6, 8, 10))+
#    coord_fixed() + 
    labs(x="PCoA 1",
         y="PCoA 2",
         color= "Vendor",
         alpha= "Day") +
    theme_classic()
}

plot_pcoa(all_pcoa) #Vegan version calculated with avgdist still looks slighltly different than mothur generated version 
plot_pcoa(exp1_pcoa)
plot_pcoa(exp2_pcoa)
plot_pcoa(`day-1_pcoa`)
plot_pcoa(`day0_pcoa`)
plot_pcoa(`day1_pcoa`)
plot_pcoa(`day2_pcoa`)
plot_pcoa(`day3_pcoa`)
plot_pcoa(`day4_pcoa`)
plot_pcoa(`day5_pcoa`)
plot_pcoa(`day6_pcoa`)
plot_pcoa(`day7_pcoa`)
plot_pcoa(`day8_pcoa`)
plot_pcoa(`day9_pcoa`)

#Check if points separate according to MiSeq run number.----
#Note, don't need to worry about 5 samples with duplicate_run entries. Since the sequences were part of the same run and were merged after make_file_fix_ids.batch with fix_ids.R
pcoa_miseq_run <- all_pcoa %>% 
  ggplot(aes(x=MDS1, y=MDS2, color = run)) +
  geom_point(size=2, alpha = 0.4) +
  coord_fixed() + 
  labs(title = NULL,
       x="PCoA 1",
       y="PCoA 2",
       color= "MiSeq Run",
       alpha= "Day") +
  theme_classic()+
  ggsave("exploratory/notebook/pcoa_by_MiSeq_run.pdf")

#Check if samples separate according to experiment number.----
#Note, don't need to worry about 5 samples with duplicate_run entries. Since the sequences were part of the same run and were merged after make_file_fix_ids.batch with fix_ids.R
pcoa_exp <- all_pcoa %>% 
  ggplot(aes(x=MDS1, y=MDS2, color = experiment)) +
  geom_point(size=2, alpha = 0.4) +
  coord_fixed() + 
  labs(title = NULL,
       x="PCoA 1",
       y="PCoA 2",
       color= "Experiment",
       alpha= "Day") +
  theme_classic()+
  ggsave("exploratory/notebook/pcoa_by_experiment.pdf")

#Statistical Analysis----
set.seed(4)
#Bray curtis distance matrix generated with vegan vegdist command using Bray-Curtis method with subsampling and iteration
all_dist <- avgdist(all, sample = 5347, distfun = vegdist, meanfun = mean, transf = NULL, iterations = 1000, dmethod = "bray")

run_adonis <- adonis(all_dist~ run, data = all_data, method = "bray")
# P = 0.001
exp_adonis <- adonis(all_dist~ experiment, data = all_pcoa, method = "bray")
# P = 0.001
vendor_adonis <- adonis(all_dist~ vendor, data = all_pcoa, method = "bray", set.seed = 4)
# P = 0.002 #Unsure why, but sometimes get 0.001 or 0.003
day_adonis <- adonis(all_dist~ day, data = all_pcoa, method = "bray")
#Df = 1, should be 10?

#With mothur amova command (based on thetayc distance matrix so values won't be an exact match):----
# Differences between experiments
# amova(phylip=data/process/vendors.subsample.thetayc.ave.dist, design=data/process/vendor.experiment.design)
# p-value < 0.001

#Differences between vendors
# amova(phylip=data/process/vendors.subsample.thetayc.ave.dist, design=data/process/vendor.vendor.design)
# P values for all the possible combinations saved as vendors.subsample.thetayc.ave.amova in test_and_old_code/amova_by_vendor

#Differences between runs
# amova(phylip=data/process/vendors.subsample.thetayc.ave.dist, design=data/process/vendor.run.design) #Note had to remove all rows with no values in the .design file to get command to work
# Fs: 2.58142 p-value: 0.006

#Differences between days
# amova(phylip=data/process/vendors.subsample.thetayc.ave.dist, design=data/process/vendor.day.design) 
# Overall p<0.001. -1 versus all other days, 0 versus all other days, 1 versus all other days, 2 vs 5-9, 3 vs 7-9
# How to interpret experiment-wise error rate and pair-wise error rate?

