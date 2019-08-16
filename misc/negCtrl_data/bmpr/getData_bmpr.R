############################################################
## Get fastq files, align, quantify and create SCE object.
## Requires:
## - sratoolkit set up and included in the path
## - R packages data.table, GEOquery, SRAdb, and Rsubread installed
## - Hisat2 installed and included in path
## - samtools installed and included in path
## - Assumes all reference-related files (indexes, GTF file) are in
## folder ./refData/. Modify as needed.
## - Gencode mm10 annotation file (release_M21)
############################################################
library(data.table)
library(GEOquery)
library(SRAdb)
library(Rsubread)

system(paste0("mkdir -p  raw"))
system(paste0("mkdir -p proc"))

############################################################
## GEO -> SRA
############################################################
gds <- getGEO("GSE70245")[[1]]
pd.dt <- data.table(pData(gds))

pd.dt <- pd.dt[grep("Bulk RNA Seq|ChIP-Seq", title, invert = T)]
pd.dt <- pd.dt[grep("poor quality", title, invert = T)] 
p1.dt <- pd.dt[characteristics_ch1 %in% c("cell type: LSK Stem cells")] 
p2.dt <- pd.dt[characteristics_ch1 %in% c("cell type: Myeloid progenitors")] 
p3.dt <- pd.dt[characteristics_ch1 %in% c("cell type: CMP Progenitors")] 
p4.dt <- pd.dt[grep("^GMP cell R", title)] 
pdRel.dt <- rbindlist(list(p1.dt, p2.dt, p3.dt, p4.dt))
dim(pdRel.dt)
info.dt <- pdRel.dt[, c("title", "geo_accession", "relation.1", "characteristics_ch1", "cell type:ch1")] 
saveRDS(info.dt, file = "info.rds")

sraIds = do.call(rbind, strsplit(as.vector(info.dt$relation.1), split = "="))[,2]
write(sraIds, file = "all_cells.txt")

############################################################
## get fastq files 
############################################################
# The fastq-dump utility from sra toolkit needs to be in path for
# below to work.  
setwd("./raw/")
cmds = paste0("nice -n19 fastq-dump --gzip --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip ",
              sraIds)
for(j in 1:length(cmds)){
    cat("Getting: ", j, "\n")
    system(cmds[j], intern = T)
}

############################################################
## align
############################################################
setwd("../")
system("./runAlignAll.sh")

############################################################
## counts
############################################################
inputBAMFiles <- Sys.glob("./proc/SRX*/hisat2/*.sorted.bam")

gtfFile <- "./refData/gencode.vM21.annotation.gtf" ## modify path as needed

counts.lst <- featureCounts(inputBAMFiles, 
                            annot.ext = gtfFile, 
                            isGTFAnnotationFile = TRUE, 
                            isPairedEnd = TRUE, 
                            nthreads = 2)
fcounts <- counts.lst$counts
colnames(fcounts) <- gsub(".sorted.bam", "", basename(inputBAMFiles))
fcounts <- fcounts[order(rownames(fcounts)), ]
saveRDS(fcounts, "./proc/all_fcounts.rds")

fannot <- counts.lst$annotation
rownames(fannot) <- fannot$GeneID
fannot <- fannot[order(rownames(fannot)), ]
saveRDS(fannot, "./proc/fc_annot.rds")

############################################################
## create procData.rds object
############################################################
library(SingleCellExperiment)
library(biomaRt)
library(stringi)

dat <- readRDS("./proc/all_fcounts.rds")
fannot <- readRDS("./proc/fc_annot.rds") 
rownames(dat) <- rownames(fannot) <- gsub("\\..*$", "", rownames(dat))
fannot.df <- data.frame(fannot$Length)
rownames(fannot.df) <- rownames(fannot)
colnames(fannot.df) <- "gene_length"

info.dt <- readRDS("./info.rds")
tt <- gsub("sra\\?term=","", basename(as.vector(info.dt$`relation.1`)))
info.dt[, "sampleId" := tt]
info.dt[, "cell_type" := gsub("cell type: ", "", characteristics_ch1)]
info.dt[, "cell_type" := gsub("Myeloid progenitors", "Mye_pro",
                              gsub("LSK Stem cells", "LSK_stem",
                                   gsub("GMP Progenitors", "GMP_pro",
                                        gsub("CMP Progenitors", "CMP_pro", cell_type))))]
info.df <- data.frame(info.dt[, c(1,6,7)])
rownames(info.df) <- info.df$sampleId
info.df <- info.df[colnames(dat), ]

mdat <- as(dat, "Matrix")
## Make single cell experiment
sce = SingleCellExperiment(assays = list(counts = mdat),
                           colData = DataFrame(info.df),
                           rowData = DataFrame(fannot.df))
sce$dbl_anno <- rep("single", ncol(sce))

sce = sce[Matrix::rowSums(counts(sce) > 0) > 0,]

#- save
saveRDS(sce, "./proc/sce_bmpr.rds")

############################################################
############################################################
############################################################
