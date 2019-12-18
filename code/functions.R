library(tidyverse)
library(cowplot)
library(magick)

### Load in metadata & create new column to give unique mouse_id based on exp. #, vendor, cage #, and mouse #.
metadata <- read_csv("data/process/vendor_metadata.v2.csv") %>% 
  unite(col = mouse_id, c(experiment, vendor, cage, mouse), remove = FALSE)
  

### Caluclate Standard Error
calc_se <- function(stdev, n){
  stdev/sqrt(n)
}

### Calculate 95% confidence intervals
lower_ci <- function(mean, se, n, conf_level = 0.95){
  lower_ci <- mean - qt(1 - ((1 - conf_level) / 2), n - 1) * se
}

upper_ci <- function(mean, se, n, conf_level = 0.95){
  upper_ci <- mean + qt(1 - ((1 - conf_level) / 2), n - 1) * se
}

#Function to have y-axis in scientific notation
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