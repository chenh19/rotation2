setwd("~/rotation2/")

if (dir.exists("./analysis/cache/")==FALSE){
  dir.create("./analysis/cache/")
}
barcode_dist_plot <- function(mtx){
  filename <- paste0("./analysis/cache/", deparse(substitute(mtx)), "_matrix_barcodes_summary.csv")
  barcode_means <- read.csv(filename, header = TRUE)
  max_barcode <- barcode_means[which.max(barcode_means$barcode_mean), "barcode_seq"]
  min_barcode <- barcode_means[which.min(barcode_means$barcode_mean), "barcode_seq"]

  # the cell with highest overall gRNAs
  target_dist <- mtx[,max_barcode]
  names(target_dist) <- 1:length(target_dist)
  filename <- paste0("./analysis/cache/1_", deparse(substitute(mtx)), "_max_barcode.tif")
  tiff(file=filename, width=8, height=5, units="in", res=600)
  barplot(target_dist, las=2, main=c(paste0("The cell with highest overall gene expression"),paste0("(barcode: ",max_barcode, ")")))
  title(xlab="Target (gene)", line=3)
  title(ylab="UMI counts", line=3.1)
  dev.off()

  # the cell with lowest overall gRNAs
  target_dist <- mtx[,min_barcode]
  names(target_dist) <- 1:length(target_dist)
  filename <- paste0("./analysis/cache/2_", deparse(substitute(mtx)), "_min_barcode.tif")
  tiff(file=filename, width=8, height=5, units="in", res=600)
  barplot(target_dist, las=2, main=c(paste0("The cell with lowest overall gene expression"),paste0("(barcode: ",min_barcode,")")))
  title(xlab="Target (gene)", line=3)
  title(ylab="UMI counts", line=3.1)
  dev.off()

  # UMI counts for all cells  (not cumulative)
  barcode_dist <- barcode_means[,"barcode_non_zero_count"]
  names(barcode_dist) <- 1:length(barcode_dist)
  filename <- paste0("./analysis/cache/3_", deparse(substitute(mtx)), "_barcode_dist.tif")
  tiff(file=filename, width=8, height=5, units="in", res=600)
  barplot(barcode_dist, las=2, main=c(paste0("Non-zero UMI counts in each barcode (cell)"),paste0("(non-cumulative)")))
  title(xlab="Barcode (Cell)", line=3.7)
  title(ylab="UMI counts", line=3.1)
  dev.off()

  # UMI counts for all cells (cdf)
  filename <- paste0("./analysis/cache/4_", deparse(substitute(mtx)), "_barcode_dist_cdf.tif")
  tiff(file=filename, width=8, height=5, units="in", res=600)
  plot(ecdf(barcode_dist), main="CDF: UMI counts in each barcode (cell)", xlab="UMI counts in each barcode (cell)", ylab="CDF")
  dev.off()

  # UMI counts for all cells (pdf)
  filename <- paste0("./analysis/cache/5_", deparse(substitute(mtx)), "_barcode_dist_pdf.tif")
  tiff(file=filename, width=8, height=5, units="in", res=600)
  PDF <- hist(barcode_dist, freq=F, breaks=150, main="PDF: UMI counts in each barcode (cell)", xlab="UMI counts in each barcode (cell)", ylab="PDF")
  dev.off()

}


Expression <- Seurat::ReadMtx(mtx = "../Morris_2021/GSM5225857_cDNA.mtx",
                                features = "../Morris_2021/GSM5225857_cDNA.features.txt",
                                cells = "../Morris_2021/GSM5225857_cDNA.barcodes.txt",
                                cell.column=1,
                                feature.column=1,
                                mtx.transpose=T)
barcode_dist_plot(Expression)
