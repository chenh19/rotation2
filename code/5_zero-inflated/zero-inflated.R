setwd("~/rotation2/")

library(Seurat)
library(dplyr)
library(ggplot2)
library(ggeasy)
library(MASS)
Expression <- readRDS("../Morris_2021/filter/Expression.rds", refhook = NULL)

if (dir.exists("./analysis/cache/")==FALSE){
  dir.create("./analysis/cache/")
}

# count means and non-zeros for the filtered Seurat subject
target_stats <- function(mtx){
  target_names <- rownames(mtx)
  target_means <- c()
  target_sums <- c()
  total <- nrow(mtx)
  i = 1
  for (target_name in target_names){
    target_means <- c(target_means, mean(mtx[target_name,]))
    target_sums <- c(target_sums, sum(mtx[target_name,] != 0 ))
    if(i%%100 == 0){
      print(paste0("Processed ",
                   format(i, big.mark=",", scientific=FALSE), " of ",
                   format(total, big.mark=",", scientific=FALSE), " targets."))
    }
    i = i+1
  }
  target_stats <- data.frame(target_names,target_means,target_sums)
  colnames(target_stats) <- c("target_name", "target_mean", "target_non_zero_count")
  filename <- paste0("./analysis/cache/Expression_filtered_means_non-zeros.csv")
  write.table(target_stats, file=filename, sep=",", row.names = FALSE)
}
target_stats(Expression@assays$RNA@counts)

# calculate zero proportions & expectations
whichmodel <- read.csv("./analysis/cache/Expression_filtered_means_non-zeros.csv", header = TRUE)
whichmodel <- dplyr::mutate(whichmodel, zero_prop = 1 - target_non_zero_count / ncol(Expression@assays$RNA@counts) )
whichmodel_10 <- dplyr::filter(whichmodel, target_mean <= 10 & zero_prop != 0)

# regression
summary(glm(zero_prop ~ target_mean, family = poisson, data = whichmodel_10)) # poisson
summary(glm.nb(zero_prop ~ target_mean, data = whichmodel_10)) # negative binomial

# plot
filename <- paste0("./analysis/cache/1_zero_proportions(Poisson).tif")
tiff(file=filename, width=8, height=5, units="in", res=600)
ggplot(data = whichmodel_10, aes(x= target_mean, y=zero_prop)) +
  geom_jitter(width=0.05, height=0.05) +
  theme_light() +
  geom_smooth(method="glm", method.args=list(family="poisson")) +
  labs(title = "Zero Proportions (Poisson regression)", x = "Gene Mean", y = "Zero Proportion") +
  ggeasy::easy_center_title()
dev.off()

filename <- paste0("./analysis/cache/2_zero_proportions(negative_binomial).tif")
tiff(file=filename, width=8, height=5, units="in", res=600)
ggplot(data = whichmodel_10, aes(x= target_mean, y=zero_prop)) +
  geom_jitter(width=0.05, height=0.05) +
  theme_light() +
  stat_smooth(method = "glm.nb", formula = y~x, color = "green") +
  labs(title = "Zero Proportions (Negative binomial regression)", x = "Gene Mean", y = "Zero Proportion") +
  ggeasy::easy_center_title()
dev.off()
