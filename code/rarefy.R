source("code/functions.R")

#Check rarefaction curves
rarefy <- read_tsv(file="data/process/vendors.subsample.rarefaction") %>% 
  select(-contains("lci-"), -contains("hci-")) %>% 
  gather(-numsampled, key=sample, value=sobs) %>% #not including numsampled, gather samples intow two columns
  mutate(sample=str_replace_all(sample, patter="0.03-", replacement="")) %>% #drop prefix 0.03
  drop_na()

#Join metadata to rarefy data frame
metadata_rarefy <- inner_join(metadata, rarefy, by = c("id" = "sample")) %>% 
  filter(vendor != "NA") 

#Plot rarefaction curves, the more parallel the curves to the x axis, the more confident you can be in the results
ggplot(metadata_rarefy, aes(x=numsampled, y=sobs, group=id, color=vendor))+
  scale_colour_manual(values=color_scheme) +
  geom_vline(xintercept=5437, color="gray", size=2) +
  geom_line()+
  coord_cartesian(xlim=c(0, 10000), ylim=c(0,200))+
  labs(subtitle="Vertical line indicates the number of sequences that samples were rarefied to",
       x="Number of sequences sampled per mouse", 
       y="Number of OTUs per mouse") +
  theme_classic()

rarefy_5437 <- metadata_rarefy %>% filter(numsampled == 5437)
#416 samples present

rarefy_5437 <- rarefy %>% filter(numsampled == 5437)