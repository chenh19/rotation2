setwd("~/rotation2/")

if (dir.exists("./analysis/cache/")==FALSE){
  dir.create("./analysis/cache/")
}

target_dist_plot <- function(mtx){
  filename <- paste0("./analysis/cache/", deparse(substitute(mtx)), "_matrix_targets_summary.csv")
  target_means <- read.csv(filename, header = TRUE)
  max_target <- target_means[which.max(target_means$target_mean), "target_name"]
  min_target <- target_means[which.min(target_means$target_mean), "target_name"]

  # the highest Hashtag in all cells
  barcode_dist <- mtx[max_target,]
  names(barcode_dist) <- 1:length(barcode_dist)
  filename <- paste0("./analysis/cache/6_", deparse(substitute(mtx)), "_max_target.tif")
  tiff(file=filename, width=8, height=5, units="in", res=600)
  barplot(barcode_dist, las=2, main=c(paste0("The highest Hashtag in all cells"),paste0("(targeting: ",max_target,")")))
  title(xlab="Barcode (Cell)", line=3.7)
  title(ylab="UMI counts", line=3.2)
  dev.off()

  # the lowest Hashtag in all cells
  barcode_dist <- mtx[min_target,]
  names(barcode_dist) <- 1:length(barcode_dist)
  filename <- paste0("./analysis/cache/7_", deparse(substitute(mtx)), "_min_target.tif")
  tiff(file=filename, width=8, height=5, units="in", res=600)
  barplot(barcode_dist, las=2, main=c(paste0("The lowest Hashtag in all cells"),paste0("(targeting: ",min_target,")")))
  title(xlab="Barcode (Cell)", line=3.7)
  title(ylab="UMI counts", line=3.2)
  dev.off()


  # UMI counts for all Hashtags  (not cumulative)
  target_dist <- target_means[,"target_non_zero_count"]
  names(target_dist) <- 1:length(target_dist)
  filename <- paste0("./analysis/cache/8_", deparse(substitute(mtx)), "_target_dist.tif")
  tiff(file=filename, width=8, height=5, units="in", res=600)
  barplot(target_dist, las=2, main=c(paste0("Non-zero UMI counts for each target (Hashtag)"),paste0("(non-cumulative)")))
  title(xlab="Target (Hashtag)", line=3)
  title(ylab="UMI counts", line=3.2)
  dev.off()

  # UMI counts for all cells (cdf)
  filename <- paste0("./analysis/cache/9_", deparse(substitute(mtx)), "_target_dist_cdf.tif")
  tiff(file=filename, width=8, height=5, units="in", res=600)
  plot(ecdf(target_dist), main="CDF: UMI counts for each target (Hashtag)", xlab="UMI counts for each target (Hashtag)", ylab="CDF")
  dev.off()

  # UMI counts for all cells (pdf)
  filename <- paste0("./analysis/cache/10_", deparse(substitute(mtx)), "_target_dist_pdf.tif")
  tiff(file=filename, width=8, height=5, units="in", res=600)
  PDF <- hist(target_dist, freq=F, breaks=100, main="PDF: UMI counts in each barcode (cell)", xlab="UMI counts in each barcode (cell)", ylab="PDF")
  dev.off()


}

Hashtag <- Seurat::ReadMtx(mtx = "../Morris_2021/GSM5225859_HTO.mtx",
                                 features = "../Morris_2021/GSM5225859_HTO.features.txt",
                                 cells = "../Morris_2021/GSM5225859_HTO.barcodes.txt",
                                 cell.column=1,
                                 feature.column=1,
                                 mtx.transpose=T)
target_dist_plot(Hashtag)
