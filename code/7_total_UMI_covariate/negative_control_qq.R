library(dplyr)

setwd("~/Developing/rotation2/")

# QQ plot
source("./code/0_prepare/qqunif.plot.R")

negative_control_naive <- read.csv("./output/SCEPTRE/negative_control_pvalues.csv")
negative_control_naive <- filter(negative_control_naive, !is.na(p_value))

negative_control_cov <- read.csv("./output/SCEPTRE/negative_control_pvalues_with_covariate.csv")
negative_control_cov <- filter(negative_control_cov, !is.na(p_value))

my.pvalue.list <- list("Negative Control Naive"= negative_control_naive$p_value,
                       "Negative Control w/ Cov"= negative_control_cov$p_value)
filename <- paste0("./analysis/cache/negative_controls_qq.tif")
tiff(file=filename, width=8, height=5, units="in", res=600)
qqunif.plot(my.pvalue.list, auto.key=list(corner=c(.95,.05)))
dev.off()
rm(list = ls())
