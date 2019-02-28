
library(Seurat)
library(SingleCellExperiment)

#- This is pretty much cut&paste from
#  https://satijalab.org/seurat/hashing_vignette.html


#---------------
#- 1. PBMC DATA
#---------------

# Load in the UMI matrix
pbmc_umi_sparse <- readRDS("./data/chpb/raw/pbmc_umi_mtx.rds")

# For generating a hashtag count matrix from fastq files, please refer to https://github.com/Hoohm/CITE-seq-Count.
# Load in the HTO count matrix
pbmc_hto <- readRDS("./data/chpb/raw/pbmc_hto_mtx.rds")

# Select cell barcodes detected by both RNA and HTO
# In the example datasets we have already filtered the cells for you, but perform this step for clarity.
joint_bcs <- intersect(colnames(pbmc_umi_sparse),colnames(pbmc_hto))

# Subset RNA and HTO counts by joint cell barcodes
pbmc_umi_sparse <- pbmc_umi_sparse[,joint_bcs]
pbmc_hto <- as.matrix(pbmc_hto[,joint_bcs])

# Confirm that the HTO have the correct names
print (rownames(pbmc_hto))

# Setup Seurat object
pbmc_hashtag <- CreateSeuratObject(raw.data = pbmc_umi_sparse)

# Normalize RNA data with log normalization
pbmc_hashtag <- NormalizeData(pbmc_hashtag,display.progress = FALSE)
# Find and scale variable genes
pbmc_hashtag <- FindVariableGenes(pbmc_hashtag,do.plot = F,display.progress = FALSE)
pbmc_hashtag <- ScaleData(pbmc_hashtag,genes.use = pbmc_hashtag@var.genes,display.progress = FALSE)

# Add HTO data as a new assay independent from RNA
pbmc_hashtag <- SetAssayData(pbmc_hashtag,assay.type = "HTO",slot = "raw.data",new.data = pbmc_hto)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
pbmc_hashtag <- NormalizeData(pbmc_hashtag,assay.type = "HTO",normalization.method = "genesCLR",display.progress = FALSE)

# If you have a very large dataset we suggest using k_function = "clara". This is a k-medoid clustering function for large applications
# You can also play with additional parameters (see documentation for HTODemux()) to adjust the threshold for classification
# Here we are using the default settings
pbmc_hashtag <- HTODemux(pbmc_hashtag,assay.type = "HTO",positive_quantile =  0.99,print.output = FALSE)
print (table(pbmc_hashtag@meta.data$hto_classification_global))

# Group cells based on the max HTO signal
pbmc_hashtag <- SetAllIdent(pbmc_hashtag,id = "hash_maxID")
pbmc_hashtag <- SetAllIdent(pbmc_hashtag,"hto_classification_global")

########################################
#- MAKE SCEs
########################################

sce.chpb = Convert(from=pbmc_hashtag, to="sce")
sce.chpb$dbl_anno_lab = sce.chpb$hto_classification_global
sce.chpb$dbl_anno = sce.chpb$hto_classification_global == "Doublet"
sce.chpb = sce.chpb[Matrix::rowSums(counts(sce.chpb)>0)>0,]

saveRDS(sce.chpb, file="./data/chpb/proc/sce_chpb.rds")
