setwd("~/Developing/rotation2/")

# Mitochondrial genes

# pull a full list of corresponding Ensembl gene IDs and symbols
#BiocManager::install("biomaRt")
library(biomaRt)
hsapiens_genes <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                        mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl"))
write.table(hsapiens_genes, file="../Morris_2021/hsapiens_genes.csv", sep=",", row.names = FALSE)
rm(hsapiens_genes)

# trim features
library(tidyr)
gene_names <- read.table("../Morris_2021/GSM5225857_cDNA.features.txt")
gene_names <- separate(gene_names, "V1", into=c("gene","version"), sep="[.]")
gene_names <- gene_names[,"gene"]
writeLines(gene_names, "../Morris_2021/cDNA_test.txt")
rm(gene_names)

# gene id to symbol conversion
library(expss)
library(dplyr)
feature <- read.table("../Morris_2021/cDNA_test.txt")
colnames(feature) <- c("id")
ref <- read.csv("../Morris_2021/hsapiens_genes.csv", na.string="")
colnames(ref) <- c("id","symbol")
ref <- mutate(ref, dup=duplicated(ref[,"id"]))
ref <- filter(ref, dup=="FALSE")
ref <- ref[,1:2]
feature <- add_columns(feature, ref, by=c("id"))
feature$gene <- ifelse(is.na(feature$symbol), feature$id, feature$symbol)
writeLines(feature$gene, "../Morris_2021/cDNA_test_feature.txt")
rm(feature,ref)

# correct for MT psudogenes
gene_card <- read.csv("../Morris_2021/GeneCards-SearchResults.csv")[,1:2] # searched "Mitochondrially Encoded" on Gene Cards and downloaded all genes
colnames(gene_card) <- c("symbol","description")
gene_card <- filter(gene_card, grepl('Mitochondrially Encoded', description) | grepl('MT-', description) | grepl('MT-', symbol))
feature <- read.csv("../Morris_2021/cDNA_test_feature.txt",header = FALSE)
colnames(feature) <- c("symbol")
feature <- add_columns(feature, gene_card, by="symbol")
feature <- mutate(feature, symbol2=(ifelse(!is.na(feature$description), paste0("MT-",feature$symbol), feature$symbol)))
feature$symbol2 <- gsub('MT-MT-','MT-',feature$symbol2)
writeLines(feature$symbol2, "../Morris_2021/cDNA_test_feature_corrected.txt")
rm(feature,gene_card)


#####################################################################################

# cDNA (filter by UMI and percent-mito)

library(Seurat)
Expression <- Seurat::ReadMtx(mtx = "../Morris_2021/GSM5225857_cDNA.mtx",
                                features = "../Morris_2021/cDNA_test_feature_corrected.txt",
                                cells = "../Morris_2021/GSM5225857_cDNA.barcodes.txt",
                                cell.column=1,
                                feature.column=1,
                                mtx.transpose=T)
Expression <- Seurat::CreateSeuratObject(counts = Expression, project = "STING-seq")
Expression[["percent.mt"]] <- Seurat::PercentageFeatureSet(Expression, pattern = "^MT-")

# plot before QC
filename <- "./analysis/cache/5_violin_before_QC.tif"
tiff(file=filename, width=8, height=5, units="in", res=1200)
Seurat::VlnPlot(Expression, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

filename <- "./analysis/cache/6_percent_mt_before_QC.tif"
tiff(file=filename, width=8, height=5, units="in", res=1200)
Seurat::FeatureScatter(Expression, feature1 = "nCount_RNA", feature2 = "percent.mt")
dev.off()

filename <- "./analysis/cache/7_nFeature_before_QC.tif"
tiff(file=filename, width=8, height=5, units="in", res=1200)
Seurat::FeatureScatter(Expression, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()


Expression <- subset(Expression, subset = nFeature_RNA >= 1400 & percent.mt < 20)
median_genes_per_cell <- median(Expression@meta.data$nFeature_RNA)
dim(Expression)
median_genes_per_cell


#####################################################################################

# HTO (filter by UMI and demuxed droplets/singlets/multiplets)

library(Seurat)
Hashtag <- Seurat::ReadMtx(mtx = "../Morris_2021/GSM5225859_HTO.mtx",
                                 features = "../Morris_2021/GSM5225859_HTO.features.txt",
                                 cells = "../Morris_2021/GSM5225859_HTO.barcodes.txt",
                                 cell.column=1,
                                 feature.column=1,
                                 mtx.transpose=T)
Hashtag <- Seurat::CreateSeuratObject(counts = Hashtag, project = "STING-seq")
Hashtag <- subset(Hashtag, subset = nFeature_RNA >= 1 & nCount_RNA < 2500)
Hashtag <- subset(x = Hashtag, cells = colnames(Hashtag)[(colnames(Hashtag) %in% c(colnames(Expression)))]) # retain only cells in Expression
dim(Hashtag)
Hashtag <- Seurat::NormalizeData(Hashtag, normalization.method = "CLR", scale.factor = 10000, margin = 2) # CLR across cells
Hashtag <- Seurat::HTODemux(Hashtag, assay = "RNA", positive.quantile = 0.91) # find singlets (I don't know how to adjust parameters for this)
Hashtag <- subset(x = Hashtag, RNA_classification.global == "Singlet")
dim(Hashtag)


#####################################################################################


# GDO (filter by UMI)

library(Seurat)
gRNA <- Seurat::ReadMtx(mtx = "../Morris_2021/GSM5225858_GDO.mtx",
                                 features = "../Morris_2021/GSM5225858_GDO.features.txt",
                                 cells = "../Morris_2021/GSM5225858_GDO.barcodes.txt",
                                 cell.column=1,
                                 feature.column=1,
                                 mtx.transpose=T)
gRNA <- Seurat::CreateSeuratObject(counts = gRNA, project = "STING-seq")
gRNA <- subset(x = gRNA, cells = colnames(gRNA)[(colnames(gRNA) %in% c(colnames(Hashtag)))]) # retain only cells in Hashtag
gRNA <- subset(gRNA, subset = nCount_RNA >= 5 & nCount_RNA <= 15000)
dim(gRNA)


#####################################################################################

# plot after QC

Expression <- subset(x = Expression, cells = colnames(Expression)[(colnames(Expression) %in% c(colnames(gRNA)))]) # retain only cells in gRNA
dim(Expression)
Hashtag <- subset(x = Hashtag, cells = colnames(Hashtag)[(colnames(Hashtag) %in% c(colnames(gRNA)))]) # retain only cells in gRNA
dim(Hashtag)
dim(gRNA)

filename <- "./analysis/cache/8_violin_after_QC.tif"
tiff(file=filename, width=8, height=5, units="in", res=1200)
Seurat::VlnPlot(Expression, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

filename <- "./analysis/cache/9_percent_mt_after_QC.tif"
tiff(file=filename, width=8, height=5, units="in", res=1200)
Seurat::FeatureScatter(Expression, feature1 = "nCount_RNA", feature2 = "percent.mt")
dev.off()

filename <- "./analysis/cache/10_nFeature_after_QC.tif"
tiff(file=filename, width=8, height=5, units="in", res=1200)
Seurat::FeatureScatter(Expression, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()


#####################################################################################

# save filtered Seurat objects as RDS for future processing

dir.create("../Morris_2021/filter", showWarnings = F)
saveRDS(Expression, "../Morris_2021/filter/Expression.rds")
saveRDS(Hashtag, "../Morris_2021/filter/Hashtag.rds")
saveRDS(gRNA, "../Morris_2021/filter/gRNAg.rds")
writeLines(gRNA@assays[["RNA"]]@counts@Dimnames[[2]], "../Morris_2021/QC_by_hang.txt")
rm(list = ls())
