setwd("~/Developing/rotation2/")
library(dplyr)
pvalue_files <- list.files(path = "./output/SCEPTRE/", pattern='*.csv')
for (pvalue_file in pvalue_files) {
  filename <- paste0("./output/SCEPTRE/", pvalue_file)
  pvalues <- read.csv(filename)
  pvalues <- filter(pvalues, !is.na(p_value))
  pvalues <- arrange(pvalues, p_value)
  print(paste0("The smallest p-values in [", pvalue_file, "] :"))
  print(head(pvalues))
}
rm(list = ls())
