library(tidyverse)

### Load in metadata
metadata <- read_csv("data/process/vendor_metadata.csv")

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
