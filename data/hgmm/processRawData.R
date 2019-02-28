library(DropletUtils)
library(Matrix)
library(readr)
library(S4Vectors)
library(biomaRt)

outputFile = "./data/hgmm/proc/sce_hgmm.rds"

#- load data sces with the hg19/mm10 counts
#==========================================
sce_hg_counts = read10xCounts("./data/hgmm/raw/hgmm_12k_raw_gene_bc_matrices_h5.h5", col.names = TRUE, group="hg19")
sce_mm_counts = read10xCounts("./data/hgmm/raw/hgmm_12k_raw_gene_bc_matrices_h5.h5", col.names = TRUE, group="mm10")

#- subset to the barcodes where we have 10X annotations hg/mm/multiplet
#======================================================================
annot           = DataFrame(read_csv("./data/hgmm/raw/analysis/gem_classification.csv"))
rownames(annot) = annot$barcode
annot           = annot[,-1]
sce_hg_counts   = sce_hg_counts[,rownames(annot)]
sce_mm_counts   = sce_mm_counts[,rownames(annot)]

colData(sce_hg_counts) = cbind(colData(sce_hg_counts),annot)
colData(sce_mm_counts) = cbind(colData(sce_mm_counts),annot)


#- map with gene symbols
hg_ens = strsplit(rowData(sce_hg_counts)$ID,split="_")
hg_ens = toupper(sapply(hg_ens, "[[",2))
mm_ens = strsplit(rowData(sce_mm_counts)$ID,split="_")
mm_ens = toupper(sapply(mm_ens, "[[",2))

## 30th Jan 2019
mart1 = useMart("ensembl", dataset="hsapiens_gene_ensembl")
mart2 = useMart("ensembl", dataset="mmusculus_gene_ensembl")

#- get the linked IDs through ensembl:
map = getLDS(attributes  = c("ensembl_gene_id"),
             filters     = "ensembl_gene_id",
             values      = hg_ens,
             mart        = mart1,
             attributesL = c("ensembl_gene_id"),
             martL       = mart2)

#- get 1:1 only
t1 = table(map[,1]); t1_uniq = names(t1)[t1==1]
t2 = table(map[,2]); t2_uniq = names(t2)[t2==1]

ind = ( map[,1] %in% t1_uniq ) & ( map[,2] %in% t2_uniq )
map = map[ind,]

map$in_hg = map[,1] %in% hg_ens
map$in_mm = map[,2] %in% mm_ens
map       = map[map$in_hg & map$in_mm,]
map$bth   = paste(map[,1],map[,2],sep="|")

hg_to_bth = map$bth; names(hg_to_bth) = map[,1]
mm_to_bth = map$bth; names(mm_to_bth) = map[,2]

hg_ind = hg_ens %in% map[,1] 
mm_ind = mm_ens %in% map[,2] 

#- ridiculous conversion orgy to end up with sparse matrix; DelayedArray does not (yet) do sparse with subset or sth....
counts_hg  = counts(sce_hg_counts); counts_hg = as(counts_hg,"matrix") ; counts_hg = Matrix(counts_hg, sparse=TRUE)
counts_mm  = counts(sce_mm_counts); counts_mm = as(counts_mm,"matrix") ; counts_mm = Matrix(counts_mm, sparse=TRUE)
counts_hg  = counts_hg[hg_ind,] #- appears actually faster to subset here only
counts_mm  = counts_mm[mm_ind,] #- appears actually faster to subset here only

#- make row and column names that have meaning
rownames(counts_hg) = hg_ens[hg_ind]; colnames(counts_hg) = colData(sce_hg_counts)$Barcode
rownames(counts_hg) = hg_to_bth[rownames(counts_hg)]
rownames(counts_mm) = mm_ens[mm_ind]; colnames(counts_mm) = colData(sce_mm_counts)$Barcode
rownames(counts_mm) = mm_to_bth[rownames(counts_mm)]
counts_mm           = counts_mm[rownames(counts_hg),]

#- add counts and stick in new sce
#=================================

counts_bth            = counts_mm + counts_hg

sce                   = sce_hg_counts
sce                   = sce[hg_ind,]
rowData(sce)$MAP      = rownames(counts_hg)
rownames(sce)         = rownames(counts_hg)
counts(sce)           = NULL
dimnames(counts_bth)  = NULL
counts(sce)           = counts_bth


call <- sce$call
call[call %in% c("hg19", "mm10")] = "Singlet"
call[call == "Multiplet"] = "Doublet"
sce$dbl_anno_lab = call
sce$dbl_anno     = sce$dbl_anno_lab == "Doublet"
sce = sce[Matrix::rowSums(counts(sce)>0)>0,]
saveRDS(sce, file = outputFile)


#- end
