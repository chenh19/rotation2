setwd("~/rotation2/")

if (dir.exists("./analysis/cache/")==FALSE){
  dir.create("./analysis/cache/")
}

target_dist_plot <- function(mtx){
  filename <- paste0("./analysis/cache/", deparse(substitute(mtx)), "_matrix_targets_summary.csv")
  target_means <- read.csv(filename, header = TRUE)
  max_target <- target_means[which.max(target_means$target_mean), "target_name"]
  min_target <- target_means[which.min(target_means$target_mean), "target_name"]

  # the highest gRNA in all cells
  barcode_dist <- mtx[max_target,]
  names(barcode_dist) <- 1:length(barcode_dist)
  filename <- paste0("./analysis/cache/6_", deparse(substitute(mtx)), "_max_target.tif")
  tiff(file=filename, width=8, height=5, units="in", res=600)
  barplot(barcode_dist, las=2, main=c(paste0("The highest expressed gene in all cells"),paste0("(gene: ",max_target,")")))
  title(xlab="Barcode (Cell)", line=3.9)
  title(ylab="UMI counts", line=3.1)
  dev.off()

  # the lowest gRNA in all cells
  barcode_dist <- mtx[min_target,]
  names(barcode_dist) <- 1:length(barcode_dist)
  filename <- paste0("./analysis/cache/7_", deparse(substitute(mtx)), "_min_target.tif")
  tiff(file=filename, width=8, height=5, units="in", res=600)
  barplot(barcode_dist, las=2, main=c(paste0("The lowest expressed gene in all cells"),paste0("(gene: ",min_target,")")))
  title(xlab="Barcode (Cell)", line=3.9)
  title(ylab="UMI counts", line=3.1)
  dev.off()


  # UMI counts for all gRNAs  (not cumulative)
  target_dist <- target_means[,"target_non_zero_count"]
  names(target_dist) <- 1:length(target_dist)
  filename <- paste0("./analysis/cache/8_", deparse(substitute(mtx)), "_target_dist.tif")
  tiff(file=filename, width=8, height=5, units="in", res=600)
  barplot(target_dist, las=2, main=c(paste0("Non-zero UMI counts for each target (gene)"),paste0("(non-cumulative)")))
  title(xlab="Target (gene)", line=3.9)
  title(ylab="UMI counts", line=3.2)
  dev.off()

  # UMI counts for all cells (cdf)
  filename <- paste0("./analysis/cache/9_", deparse(substitute(mtx)), "_target_dist_cdf.tif")
  tiff(file=filename, width=8, height=5, units="in", res=600)
  plot(ecdf(target_dist), main="CDF: UMI counts for each target (gene)", xlab="UMI counts for each target (gene)", ylab="CDF")
  dev.off()

  # UMI counts for all cells (pdf)
  filename <- paste0("./analysis/cache/10_", deparse(substitute(mtx)), "_target_dist_pdf.tif")
  tiff(file=filename, width=8, height=5, units="in", res=600)
  PDF <- hist(target_dist, freq=F, breaks=150, main="PDF: UMI counts in each barcode (cell)", xlab="UMI counts in each barcode (cell)", ylab="PDF")
  dev.off()


}

Expression <- Seurat::ReadMtx(mtx = "../Morris_2021/GSM5225857_cDNA.mtx",
                                features = "../Morris_2021/GSM5225857_cDNA.features.txt",
                                cells = "../Morris_2021/GSM5225857_cDNA.barcodes.txt",
                                cell.column=1,
                                feature.column=1,
                                mtx.transpose=T)
target_dist_plot(Expression)
