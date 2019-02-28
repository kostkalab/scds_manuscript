
library(Seurat)
library(SingleCellExperiment)
#- This is pretty much cut&paste from
#  https://satijalab.org/seurat/hashing_vignette.html

#-------------------
#- 1. CELL LINE DATA
#-------------------

# Read in UMI count matrix for RNA
hto12_umi <- readRDS("./data/chcl/raw/hto12_umi_mtx.rds")

# Read in HTO count matrix
hto12_hto_raw <- readRDS("./data/chcl/raw/hto12_hto_mtx.rds")

# Select cell barcodes detected in both RNA and HTO
cells_use <- intersect(rownames(hto12_hto_raw),colnames(hto12_umi))

# Create Seurat object and add HTO data
hto12 <- CreateSeuratObject(hto12_umi[,cells_use],min.genes = 300)
hto12 <- SetAssayData(hto12,assay.type = "HTO",slot = "raw.data",new.data = t(hto12_hto_raw[hto12@cell.names,1:12]))

# Normalize data
hto12 <- NormalizeData(hto12,display.progress = FALSE)
hto12 <- NormalizeData(hto12,assay.type = "HTO",normalization.method = "genesCLR",display.progress = FALSE)

#- Demultiplex
hto12 <- HTODemux(hto12,assay.type = "HTO",positive_quantile =  0.99, print.output = TRUE)


########################################
#- MAKE SCEs
########################################

sce.chcl              = Convert(from=hto12, to="sce")
sce.chcl$dbl_anno_lab = sce.chcl$hto_classification_global
sce.chcl$dbl_anno     = sce.chcl$hto_classification_global == "Doublet"
sce.chcl              = sce.chcl[Matrix::rowSums(counts(sce.chcl)>0)>0,]

saveRDS(sce.chcl, file="./data/chcl/proc/sce_chcl.rds")
