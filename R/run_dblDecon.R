library(Seurat)
library(DoubletDecon)
library(dplyr)
library(SingleCellExperiment)

source("./R/dd/fixed_dd.R")

run_dblDecon <- function(sce, species = "hsa"){

  if(!dir.exists("./results/tmp/dd")) dir.create("./results/tmp/dd", recursive = TRUE)
  location = "./results/tmp/dd/"
  filename = "dd_out"

  ## create fully-processed Seurat object
  seu <- CreateSeuratObject(assay(sce, "counts"))
  seu <- NormalizeData(seu)
  seu <- ScaleData(seu, vars.to.regress = "nUMI") ## similar to processing for DoubletFinder input
  seu <- FindVariableGenes(seu, do.plot = FALSE)

  seu <- RunPCA(seu, pc.genes = seu@var.genes, pcs.print = 0)
  seu <- FindClusters(seu)
  seu.markers <- FindAllMarkers(object = seu)
  
  seu.markers = seu.markers[!grepl("\\.",rownames(seu.markers)),]
  seu.markers <- seu.markers %>% group_by(cluster) %>% top_n(50, avg_logFC)

  ## Seurat output files for pre-processing with their function
  expression <- as(attributes(seu)$data, "matrix")
  clusters <- seu@ident

  exprsFile <- paste0(location, "input_exprs.txt")
  clustersFile <- paste0(location, "input_clusters.txt")
  genesFile <- paste0(location, "input_genes.txt")
  write.table(expression, file = exprsFile, sep = "\t", row.names=T, col.names = NA, quote=F)
  write.table(clusters, file = clustersFile, sep = "\t", row.names=T, col.names = NA, quote=F)
  write.table(seu.markers, file = genesFile, sep = "\t", row.names=T, col.names = NA, quote=F)

  newFiles <- Seurat_Pre_Process(exprsFile, genesFile, clustersFile)

  ## run DoubletDecon
  pdf("/dev/null")
  rawDataFile <- newFiles$newExpressionFile
  groupsFile <- newFiles$newGroupsFile
  results = myMain_Doublet_Decon(rawDataFile = rawDataFile,
                               groupsFile = groupsFile,
                               filename = filename,
                               location = location,
                               removeCC = FALSE,
                               species = species,
                               rhop = 1,
                               write = TRUE,
                               heatmap = FALSE,
                               centroids = TRUE,
                               downsample = "none",
                               sample_num = NULL)
  dev.off()
  saveRDS(results, file = paste0(location, "results.rds"))
  doublets_res <- results$Final_doublets_groups
  colData(sce)$dblDecon_call = rep("Singlet", ncol(sce))
  colData(sce)[rownames(doublets_res), ]$dblDecon_call = "Doublet"
  metadata(sce)$dbl_decon_results = results
  return(sce)
}