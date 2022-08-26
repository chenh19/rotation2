setwd("~/Developing/rotation2/")
library(tidyverse)
library(cowplot)
library(Matrix)
library(sceptre)
side <- "left"

# read gene and gRNA matrices
Expression <- readRDS("../Morris_2021/filter/Expression.rds", refhook = NULL)
gRNA <- readRDS("../Morris_2021/filter/gRNA.rds", refhook = NULL)
gene_matrix2 <- as(Expression@assays$RNA@counts, "dgTMatrix")
gRNA_matrix2 <- as(gRNA@assays$RNA@counts, "dgTMatrix")
rm("Expression","gRNA")

# avoid Matrix::t() dimension issue
gRNA_null <- Matrix(nrow = 12, ncol = ncol(gRNA_matrix2), data = 0, sparse = TRUE)
gRNA_null <- as(gRNA_null, "dgTMatrix")
rownames(gRNA_null) <- c("NTC-1-null","NTC-2-null","NTC-3-null","NTC-4-null","NTC-5-null","NTC-6-null","NTC-7-null","NTC-8-null","NTC-9-null","NTC-10-null","NTC-11-null","NTC-12-null")
colnames(gRNA_null) <- colnames(gRNA_matrix2)
gRNA_matrix2 <- rbind(gRNA_matrix2, gRNA_null)
gRNA_matrix2 <- as(gRNA_matrix2, "dgTMatrix")

# generate covariate matrix
covariate_matrix2 <- as.data.frame(colnames(gene_matrix2))
rownames(covariate_matrix2) <- colnames(gene_matrix2)
covariate_matrix2$lg_gRNA_lib_size <- rnorm(nrow(covariate_matrix2), 10, 0.5)
covariate_matrix2$lg_gene_lib_size <- rnorm(nrow(covariate_matrix2), 10, 0.3)
covariate_matrix2$p_mito <- rnorm(nrow(covariate_matrix2), 0.01, 0.001)
covariate_matrix2$batch <- as.factor(ifelse(covariate_matrix2$lg_gRNA_lib_size >= 10, paste0("prep_batch_1"), paste0("prep_batch_2")))
covariate_matrix2 <- covariate_matrix2[c("lg_gRNA_lib_size","lg_gene_lib_size",
                                         "p_mito","batch")]

# generate perturbation matrix
perturbation_matrix2 <- threshold_gRNA_matrix(gRNA_matrix2)

# subset the expression matrix by gene
genes <- lapply(paste0("./output/SCEPTRE/gene-gRNA_pairs/",
                list.files(path = "./output/SCEPTRE/gene-gRNA_pairs/", pattern='*_pairs.csv')),
                read.csv)
genes <- do.call(rbind.data.frame, genes)[,"Gene"]
Sys.sleep(3)
gene_matrix2 <- gene_matrix2[(rownames(gene_matrix2) %in% genes), ]
rm("genes")

# gRNA group table
gRNA_groups_table2 <- as_tibble(read.csv("./output/SCEPTRE/gene-gRNA_pairs/gRNA_groups_table2.csv", header = TRUE))

# combine perturbation
combined_perturbation_matrix2 <- combine_perturbations(perturbation_matrix = perturbation_matrix2,
                                                      gRNA_groups_table = gRNA_groups_table2)

# gene gRNA group pairs
gene_gRNA_group_pairs2 <- as.data.frame(read.csv("./output/SCEPTRE/gene-gRNA_pairs/gene_gRNA_group_pairs2.csv", header = TRUE)[,1:3])
gene_gRNA_group_pairs2$check <- paste0(gene_gRNA_group_pairs2$gene_id, gene_gRNA_group_pairs2$gRNA_group, gene_gRNA_group_pairs2$pair_type)
gene_gRNA_group_pairs2 <- gene_gRNA_group_pairs2[!duplicated(gene_gRNA_group_pairs2[,"check"]),]
gene_gRNA_group_pairs2 <- gene_gRNA_group_pairs2[gene_gRNA_group_pairs2$gene_id %in% rownames(gene_matrix2),]
gene_gRNA_group_pairs2 <- as_tibble((gene_gRNA_group_pairs2)[,1:3])

# final compute
gene_gRNA_group_pairs2$pair_type <- as.factor(gene_gRNA_group_pairs2$pair_type)
colnames(gene_gRNA_group_pairs2) <- colnames(gene_gRNA_group_pairs)
gene_matrix <- gene_matrix2
gRNA_matrix <- gRNA_matrix2
combined_perturbation_matrix <- combined_perturbation_matrix2
covariate_matrix <- covariate_matrix2
gene_gRNA_group_pairs <- gene_gRNA_group_pairs2

result2 <- run_sceptre_high_moi(gene_matrix = gene_matrix,
                               combined_perturbation_matrix = combined_perturbation_matrix,
                               covariate_matrix = covariate_matrix,
                               gene_gRNA_group_pairs = gene_gRNA_group_pairs,
                               side = side)
write.table(result2, file="./analysis/cache/result2.csv", sep=",", row.names = FALSE)
