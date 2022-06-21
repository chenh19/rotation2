setwd("~/rotation2/")

if (dir.exists("./analysis/cache/")==FALSE){
  dir.create("./analysis/cache/")
}

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
  filename <- paste0("./analysis/cache/", deparse(substitute(mtx)), "_matrix_targets_summary.csv")
  write.table(target_stats, file=filename, sep=",", row.names = FALSE)
}


Expression <- Seurat::ReadMtx(mtx = "../Morris_2021/GSM5225857_cDNA.mtx",
                                features = "../Morris_2021/GSM5225857_cDNA.features.txt",
                                cells = "../Morris_2021/GSM5225857_cDNA.barcodes.txt",
                                cell.column=1,
                                feature.column=1,
                                mtx.transpose=T)
target_stats(Expression)
