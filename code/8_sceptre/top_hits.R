# top hits
setwd("~/Developing/rotation2/")
library(dplyr)
sceptre <- read.csv("./output/SCEPTRE/result2.csv")
hits <- filter(sceptre, !is.na(p_value) & gRNA_type == "candidate" & p_value < 0.05/nrow(sceptre))
write.table(hits, file="./analysis/cache/hits.csv", sep=",", row.names = FALSE)

# cleanup
rm(list = ls())
