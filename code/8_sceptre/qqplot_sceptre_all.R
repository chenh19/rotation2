setwd("~/Developing/rotation2/")
library(dplyr)

# QQ plot negative controls
source("./code/0_prepare/qqunif.plot.R")
negative_control_naive <- read.csv("./output/SCEPTRE/negative_control_pvalues.csv")
negative_control_naive <- filter(negative_control_naive, !is.na(p_value))
negative_control_cov <- read.csv("./output/SCEPTRE/negative_control_pvalues_with_covariate.csv")
negative_control_cov <- filter(negative_control_cov, !is.na(p_value))
negative_control_sceptre <- read.csv("./output/SCEPTRE/result2.csv")
negative_control_sceptre <- filter(negative_control_sceptre, !is.na(p_value))
negative_control_sceptre <- filter(negative_control_sceptre, gRNA_type == "negative_control")

my.pvalue.list <- list("Negative Control (Naive NB)"= negative_control_naive$p_value,
                       "Negative Control (NB w/ Cov)"= negative_control_cov$p_value,
                       "Negative Control (SCEPTRE)"= negative_control_sceptre$p_value)
filename <- paste0("./analysis/cache/negative_controls_qq_sceptre.tif")
tiff(file=filename, width=8, height=5, units="in", res=600)
qqunif.plot(my.pvalue.list, auto.key=list(corner=c(.95,.05)))
dev.off()


# QQ plot sceptre
source("./code/0_prepare/qqunif.plot.R")
sceptre <- read.csv("./output/SCEPTRE/result2.csv")
negative_control_sceptre <- filter(sceptre, !is.na(p_value) & gRNA_type == "negative_control")
positive_control_sceptre <- filter(sceptre, !is.na(p_value) & gRNA_type == "positive_control")
candidate_sceptre <- filter(sceptre, !is.na(p_value) & gRNA_type == "candidate")

my.pvalue.list <- list("Negative Control (SCEPTRE)"= negative_control_sceptre$p_value,
                       "Positive Control (SCEPTRE)"= positive_control_sceptre$p_value,
                       "Candidate (SCEPTRE)"= candidate_sceptre$p_value)
filename <- paste0("./analysis/cache/qq_sceptre_all.tif")
tiff(file=filename, width=8, height=5, units="in", res=600)
qqunif.plot(my.pvalue.list, auto.key=list(corner=c(.95,.05)))
dev.off()

# cleanup
rm(list = ls())
