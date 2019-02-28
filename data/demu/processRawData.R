library(BiocFileCache)
library(data.table)
library(SingleCellExperiment)


bfc <- BiocFileCache("./data/demu/raw", ask = FALSE)

base.path <- "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2560nnn/GSM2560248/suppl"
barcode.fname <- bfcrpath(bfc, file.path(base.path,
                                         "GSM2560248%5Fbarcodes%2Etsv%2Egz"))
counts.fname <- bfcrpath(bfc, file.path(base.path,
                                        "GSM2560248%5F2%2E1%2Emtx%2Egz"))
## downloaded separately from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96583
## GSM2560248
gene.fname <- bfcrpath(bfc, file.path("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE96nnn/GSE96583/suppl/GSE96583%5Fbatch2%2Egenes%2Etsv%2Egz"))

library(scater)
library(Matrix)
gene.info <- read.table(gene.fname, stringsAsFactors=FALSE)
colnames(gene.info) <- c("Ensembl", "Symbol")
sce <- SingleCellExperiment(
  list(counts = as(readMM(counts.fname), "dgCMatrix")),
  rowData = gene.info,
  colData = DataFrame(Barcode = readLines(barcode.fname))
)
rownames(sce) <- uniquifyFeatureNames(
  rowData(sce)$Ensembl, rowData(sce)$Symbol)
colnames(sce) <- sce$Barcode
sce

## read in annotation
annot = DataFrame(fread("./data/demu/raw/ye1.ctrl.8.10.sm.best"))
rownames(annot) = annot$BARCODE
annot = annot[colnames(sce),]
colData(sce) = cbind(colData(sce),annot)

## get doublet true labels
best <- sce$BEST
best <- substr(sce$BEST, start = 1, stop = 3)
best[best == "DBL"] = "Doublet"
best[best == "SNG"] = "Singlet"
best[best == "AMB"] = "Ambiguous"

sce$dbl_anno_lab = best
sce$dbl_anno     = sce$dbl_anno_lab == "Doublet"
sce = sce[Matrix::rowSums(counts(sce)>0)>0,]

saveRDS(sce, "./data/demu/proc/sce_demu.rds")
