#install.packages("devtools")
#devtools::install_github("katsevich-lab/sceptre")
#install.packages("tidyverse")
#install.packages("cowplot")

# SCEPTRE example
library(tidyverse)
library(cowplot)
library(Matrix)
library(sceptre)

data(gene_matrix)
data(gRNA_matrix)
data(covariate_matrix)
perturbation_matrix <- threshold_gRNA_matrix(gRNA_matrix)
data("gRNA_groups_table")
combined_perturbation_matrix <- combine_perturbations(perturbation_matrix = perturbation_matrix,
                                                      gRNA_groups_table = gRNA_groups_table)

data(gene_gRNA_group_pairs)
side <- "left"
result <- run_sceptre_high_moi(gene_matrix = gene_matrix,
                               combined_perturbation_matrix = combined_perturbation_matrix,
                               covariate_matrix = covariate_matrix,
                               gene_gRNA_group_pairs = gene_gRNA_group_pairs,
                               side = side)
