setwd("~/Developing/rotation2/")

if (dir.exists("./analysis/cache/")==FALSE){
  dir.create("./analysis/cache/")
}
barcode_dist_plot_before_QC <- function(mtx){

  filename <- paste0("./analysis/cache/", deparse(substitute(mtx)), "_matrix_barcodes_summary.csv")
  barcode_means <- read.csv(filename, header = TRUE)
  barcode_means <- dplyr::filter(barcode_means, barcode_non_zero_count >=50) # & barcode_non_zero_count <=4000


  # UMI counts for all cells  (not cumulative)
  barcode_dist <- barcode_means[,"barcode_non_zero_count"]
  names(barcode_dist) <- 1:length(barcode_dist)

  # UMI counts for all cells (cdf)
  filename <- paste0("./analysis/cache/1_", deparse(substitute(mtx)), "_barcode_dist_qc_cdf.tif")
  tiff(file=filename, width=8, height=5, units="in", res=600)
  plot(ecdf(barcode_dist), main="CDF: UMI counts in each barcode (cell)", xlab="UMI counts in each barcode (cell)", ylab="CDF")
  abline(v=850, col="darkgray", lwd=2, lty=2)
  abline(v=1400, col="red", lwd=2, lty=2)
  text(1900, 0.5, "New QC: UMI count >= 1400", col="red", adj=0)
  text(1900, 0.4, "Original QC: UMI count >= 850", col="darkgray", adj=0)
  text(1900, 0.3, "(UMI count < 50 not plotted)", col="lightgray", adj=0)
  dev.off()

  # UMI counts for all cells (pdf)
  filename <- paste0("./analysis/cache/2_", deparse(substitute(mtx)), "_barcode_dist_qc_pdf.tif")
  tiff(file=filename, width=8, height=5, units="in", res=600)
  PDF <- hist(barcode_dist, freq=F, breaks=150, main="PDF: UMI counts in each barcode (cell)", xlab="UMI counts in each barcode (cell)", ylab="PDF")
  abline(v=850, col="darkgray", lwd=2, lty=2)
  abline(v=1400, col="red", lwd=2, lty=2)
  text(1900, 0.003, "New QC: UMI count >= 1400", col="red", adj=0)
  text(1900, 0.0025, "Original QC: UMI count >= 850", col="darkgray", adj=0)
  text(1900, 0.002, "(UMI count < 50 not plotted)", col="lightgray", adj=0)
  dev.off()

}

barcode_dist_plot_after_QC <- function(mtx){

  filename <- paste0("./analysis/cache/", deparse(substitute(mtx)), "_matrix_barcodes_summary.csv")
  barcode_means <- read.csv(filename, header = TRUE)
  barcode_means <- dplyr::filter(barcode_means, barcode_non_zero_count >=1400)


  # UMI counts for all cells  (not cumulative)
  barcode_dist <- barcode_means[,"barcode_non_zero_count"]
  names(barcode_dist) <- 1:length(barcode_dist)

  # UMI counts for all cells (cdf)
  filename <- paste0("./analysis/cache/3_", deparse(substitute(mtx)), "_barcode_dist_qc_cdf.tif")
  tiff(file=filename, width=8, height=5, units="in", res=600)
  plot(ecdf(barcode_dist), main="CDF: UMI counts in each barcode (cell)", xlab="UMI counts in each barcode (cell)", ylab="CDF")
  #abline(v=850, col="red", lwd=2, lty=2)
  text(5500, 0.6, "(After new QC: UMI counts >= 1400)", col="blue", adj=0)
  #text(200, 0.9, "(UMI count < 50 excluded)")
  dev.off()

  # UMI counts for all cells (pdf)
  filename <- paste0("./analysis/cache/4_", deparse(substitute(mtx)), "_barcode_dist_qc_pdf.tif")
  tiff(file=filename, width=8, height=5, units="in", res=600)
  PDF <- hist(barcode_dist, freq=F, breaks=150, main="PDF: UMI counts in each barcode (cell)", xlab="UMI counts in each barcode (cell)", ylab="PDF")
  #abline(v=850, col="red", lwd=2, lty=2)
  text(5500, 3e-4, "(After new QC: UMI counts >= 1400)", col="blue", adj=0)
  #text(290, 0.006, "(UMI count < 50 excluded)")
  dev.off()

}

barcode_dist_plot_before_QC(Expression)
barcode_dist_plot_after_QC(Expression)
