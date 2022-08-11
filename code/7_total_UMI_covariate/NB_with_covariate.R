setwd("~/Developing/rotation2/")

library (MASS)
library(Matrix)
library(dplyr)

# try function
try <- function(code, silent = FALSE) {
    tryCatch(code, error = function(c) {
    if (!silent) {"Error Message"}
    else{code}})}


# load data
Expression <- readRDS("../Morris_2021/filter/Expression.rds", refhook = NULL)
gRNA <- readRDS("../Morris_2021/filter/gRNA.rds", refhook = NULL)
gene_matrix2 <- Expression@assays$RNA@counts
gRNA_matrix2 <- gRNA@assays$RNA@counts
rm("Expression","gRNA")

# trim genes
genes <- lapply(paste0("./output/SCEPTRE/gene-gRNA_pairs/",
                list.files(path = "./output/SCEPTRE/gene-gRNA_pairs/", pattern='*_pairs.csv')),
                read.csv)
genes <- do.call(rbind.data.frame, genes)[,"Gene"]
Sys.sleep(3)
gene_matrix2 <- gene_matrix2[(rownames(gene_matrix2) %in% genes), ]
rm("genes")


# NB regression
nbregression <- function(sample){
  filename <- paste0("./output/SCEPTRE/gene-gRNA_pairs/", deparse(substitute(sample)), "_pairs.csv")
  pairs <- read.csv(filename)
  pairs <- pairs[(pairs$Gene %in% rownames(gene_matrix2)),]

  genes <- c()
  gRNAs <- c()

  count <- nrow(pairs)
  pvalues <- c()
  i = 1
  for (i in 1:count){
    gene <- pairs[i,1]
    gene_data <- gene_matrix2[gene,]
    gRNA <- pairs[i,2]
    gRNA_data <- gRNA_matrix2[gRNA,]
    total_UMI <- colSums(gene_matrix2)
    pair <- data.frame(gene_data, gRNA_data, total_UMI)
    test <- try(glm.nb(gene_data ~ gRNA_data + total_UMI, data = pair, start = NULL))
    genes <- c(genes, gene)
    gRNAs <- c(gRNAs, gRNA)

    if (is.list(test)){
      pvalue <- coef(summary(test))[2,4] # negative binomial
      pvalues <- c(pvalues, pvalue)
    } else {
      pvalues <- c(pvalues, NA)
    }

    if (count >= 100) {
      if (i%%100 == 0) {
          print(paste0("Processed ",
                 format(i, big.mark=",", scientific=FALSE), " of ",
                 format(count, big.mark=",", scientific=FALSE), " ",
                 "[", deparse(substitute(sample)), "] pairs."))
      }
    }

    i = i+1
  }
  output <- data.frame(genes, gRNAs, pvalues)
  colnames(output) <- c("Gene", "gRNA", "p_value")
  filename <- paste0("./output/SCEPTRE/", deparse(substitute(sample)), "_pvalues_with_covariate.csv")
  write.table(output, file=filename, sep=",", row.names = FALSE)
  print(paste0("Processed all ",format(count, big.mark=",", scientific=FALSE), " ",
                 "[", deparse(substitute(sample)), "] pairs."))
}

nbregression(positive_control)
nbregression(candidate)
nbregression(negative_control)


# QQ plot
source("./code/0_prepare/qqunif.plot.R")
negative_control <- read.csv("./output/SCEPTRE/negative_control_pvalues_with_covariate.csv")
negative_control <- filter(negative_control, !is.na(p_value))
positive_control <- read.csv("./output/SCEPTRE/positive_control_pvalues_with_covariate.csv")
positive_control <- filter(positive_control, !is.na(p_value))
candidate <- read.csv("./output/SCEPTRE/candidate_pvalues_with_covariate.csv")
candidate <- filter(candidate, !is.na(p_value))
my.pvalue.list <- list("Negative Control"= negative_control$p_value,
                       "Positive Control"= positive_control$p_value,
                       "Candidate"= candidate$p_value)
filename <- paste0("./analysis/cache/NB_with_covarite_qq.tif")
tiff(file=filename, width=8, height=5, units="in", res=600)
qqunif.plot(my.pvalue.list, auto.key=list(corner=c(.95,.05)))
dev.off()

# cleanup
rm(list = ls())
