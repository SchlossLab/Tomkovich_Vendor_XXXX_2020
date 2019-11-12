source("code/functions.R")

pcoa_data <- read_tsv("data/process/vendors.subsample.thetayc.ave.pcoa.axes")

pcoa_plot <- right_join(metadata, pcoa_data, by=c("id" = "group")) %>%
  ggplot(aes(x=axis1, y=axis2, color=vendor)) +
	geom_point(shape=19, size=2) +
	scale_color_manual(name=NULL,
		values=c("blue", "red", "black", "green", "purple", "orange"),
		breaks=c("Jackson", "Charles River", "Taconic", 
		         "Schloss", "Young", "Envigo"),
		labels=c("Jackson Laboratories", "Charles River Labs", "Taconic Biosciences",
		         "Schloss Colony", "Young Colony", "Envigo")) +
	coord_fixed() + 
	labs(x="PCo Axis 1", y="PCo Axis 2") +
	theme_classic()

pcoa_plot
