setwd("~/rotation2/")

if (dir.exists("./analysis/cache/")==FALSE){
  dir.create("./analysis/cache/")
}

barcode_stats <-function(mtx){
  barcode_seqs <- colnames(mtx)
  barcode_means <- c()
  barcode_sums <- c()
  total <- ncol(mtx)
  i = 1
  for (barcode_seq in barcode_seqs){
    barcode_means <- c(barcode_means, mean(mtx[,barcode_seq]))
    barcode_sums <- c(barcode_sums, sum(mtx[,barcode_seq] != 0 ))
    if(i%%100 == 0){
    print(paste0("Processed ",
                 format(i, big.mark=",", scientific=FALSE), " of ",
                 format(total, big.mark=",", scientific=FALSE), " barcodes."))
    }
    i = i+1
  }
  barcode_stats <- data.frame(barcode_seqs,barcode_means,barcode_sums)
  colnames(barcode_stats) <- c("barcode_seq", "barcode_mean", "barcode_non_zero_count")
  filename <- paste0("./analysis/cache/", deparse(substitute(mtx)), "_matrix_barcodes_summary.csv")
  write.table(barcode_stats, file=filename, sep=",", row.names = FALSE)
}


Expression <- Seurat::ReadMtx(mtx = "../Morris_2021/GSM5225857_cDNA.mtx",
                                features = "../Morris_2021/GSM5225857_cDNA.features.txt",
                                cells = "../Morris_2021/GSM5225857_cDNA.barcodes.txt",
                                cell.column=1,
                                feature.column=1,
                                mtx.transpose=T)
barcode_stats(Expression)
