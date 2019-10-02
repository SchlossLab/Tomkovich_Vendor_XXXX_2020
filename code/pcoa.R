library(tidyverse)
library(readxl)

pcoa_data <- read_tsv(file="data/mothur/.pcoa.axes")
metadata <- 
meadata <- inner_join(metadata, pcoa, by=c('sample'='group'))

ggplot(metadata_pcoa, aes(x=axis1, y=axis2, color=vendor)) +
	geom_point(shape=19, size=2) +
	scale_color_manual(name=NULL,
		values=c("blue", "red", "black", "green", "purple", "orange"),
		breaks=c("jackson", "charles_river", "tac"),
		labels=c("Jackson Laboratories",) +
	coord_fixed() + 
	labs(x="PCo Axis 1", y="PCo Axis 2") +
	theme_classic()


