############################################################
############################################################
library(data.table)
library(SingleCellExperiment)
library(biomaRt)
library(stringi)

dat <- fread("./raw/GSE100597_count_table_QC_filtered.txt.gz", head = T, stringsAsFactors = F)
gns <- gns_dat <- dat$V1

dat <- dat[, -1]
colnames(dat) <- gsub("Single", "single", gsub("Double", "double", gsub("Cell_", "Cell", colnames(dat))))

## excluding VE cells
cns <- colnames(dat)
cns <- cns[!grepl("VE", cns)]
dat <- dat[, ..cns]

col_ann <- unlist(lapply(strsplit(colnames(dat), split = "_"), function(x){
    return(x[length(x)])
}))

mdat <- as(dat, "Matrix")
colnames(mdat) <- cns
rownames(mdat) <- gns

## Make single cell experiment
sce = SingleCellExperiment(assays = list(counts = mdat),
                           colData = DataFrame(col_ann),
                           rowData = DataFrame(gns))
sce$dbl_anno <- as.vector(sce$col_ann)
sce = sce[, sce$dbl_anno != "double"] ## keeping only singlets

sce = sce[Matrix::rowSums(counts(sce) > 0) > 0,]

#- save
saveRDS(sce, "./proc/sce_micm.rds")

############################################################
############################################################
############################################################




