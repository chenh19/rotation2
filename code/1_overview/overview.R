setwd("~/rotation2/")

overall <- function(mtx){
  cat(paste0('The [', deparse(substitute(mtx)), '] matrix has: \n',
        '- ', format(mtx@Dim[1], big.mark=",", scientific=FALSE),
        ' rows/genes/targets \n',
        '- ', format(mtx@Dim[2], big.mark=",", scientific=FALSE),
        ' columns/barcodes/cells \n',
        '- ', format(length(mtx), big.mark=",", scientific=FALSE),
        ' values in total \n',
        '- ', format(sum(mtx != 0), big.mark=",", scientific=FALSE),
        ' values that are non-zero \n',
        '- ', format(sum(mtx == 1), big.mark=",", scientific=FALSE),
        ' values that are 1 \n',
        '- ', format(sum(mtx > 1), big.mark=",", scientific=FALSE),
        ' values that are bigger than 1 \n',
        '- ', format(sum(mtx > 10), big.mark=",", scientific=FALSE),
        ' values that are bigger than 10 \n',
        '- ', format(sum(mtx > 100), big.mark=",", scientific=FALSE),
        ' values that are bigger than 100 \n',
        '- ', format(sum(mtx > 1000), big.mark=",", scientific=FALSE),
        ' values that are bigger than 1,000 \n',
        '- ', format(sum(mtx > 10000), big.mark=",", scientific=FALSE),
        ' values that are bigger than 10,000 \n',
        '- ', format(sum(mtx > 100000), big.mark=",", scientific=FALSE),
        ' values that are bigger than 100,000 \n',
        ' \n'))
}

# load expression data
Expression <- Seurat::ReadMtx(mtx = "../Morris_2021/GSM5225857_cDNA.mtx",
                                features = "../Morris_2021/GSM5225857_cDNA.features.txt",
                                cells = "../Morris_2021/GSM5225857_cDNA.barcodes.txt",
                                cell.column=1,
                                feature.column=1,
                                mtx.transpose=T)

# load gRNA data
gRNA <- Seurat::ReadMtx(mtx = "../Morris_2021/GSM5225858_GDO.mtx",
                                 features = "../Morris_2021/GSM5225858_GDO.features.txt",
                                 cells = "../Morris_2021/GSM5225858_GDO.barcodes.txt",
                                 cell.column=1,
                                 feature.column=1,
                                 mtx.transpose=T)

# load Hashtag data
Hashtag <- Seurat::ReadMtx(mtx = "../Morris_2021/GSM5225859_HTO.mtx",
                                 features = "../Morris_2021/GSM5225859_HTO.features.txt",
                                 cells = "../Morris_2021/GSM5225859_HTO.barcodes.txt",
                                 cell.column=1,
                                 feature.column=1,
                                 mtx.transpose=T)

overall(Expression)
overall(gRNA)
overall(Hashtag)
