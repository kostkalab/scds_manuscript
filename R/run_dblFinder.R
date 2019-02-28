library(scran)
library(Seurat)
library(DoubletFinder)

run_dblFinder <- function(sce, use_celltype_anno = FALSE, gitParamSettings = TRUE){
  pN = 0.25

  ## create fully-processed Seurat object
  seu <- CreateSeuratObject(assay(sce, "counts"))
  seu <- NormalizeData(seu)
  seu <- ScaleData(seu, vars.to.regress = "nUMI")

  ## if use the cutoffs as in the github example code
  if(gitParamSettings){
    seu <- FindVariableGenes(seu, x.low.cutoff = 0.0125, y.cutoff = 0.25, do.plot = FALSE)
  } else {
    seu <- FindVariableGenes(seu, do.plot = FALSE)
  }

  seu <- RunPCA(seu, pc.genes = seu@var.genes, pcs.print = 0)
  seu <- RunTSNE(seu, dims.use = 1:10, verbose=TRUE)
  ## ASB added:
  seu <- FindClusters(seu, reduction.type = "pca", dims.use = 1:10) ## , resolution = 0.6, print.output = 0, save.SNN = TRUE)
  ## without dims.use set, the paramSweep does not work

  ## pK Identification
  sweep.res.list <- paramSweep(seu)
  sweep.stats <- summarizeSweep(sweep.res.list) ## GT = FALSE is default

  #- can't switch off plotting
  pdf("/dev/null")
  bcmvn <- find.pK(sweep.stats)
  dev.off()
  ## choose pK per best practices outlined on github

  pK <- as.numeric(as.vector(bcmvn[which.max(bcmvn$BCmetric), ]$pK))

## expected number of doublets per Poisson assumption
  nExp_poi <- round(0.075 * length(seu@cell.names))
  cn <- paste(pN, pK, nExp_poi, sep = "_")
  cat(cn, "\n")

  ## run DoubletFinder without prior cluster annotation information - default
  seu <- doubletFinder(seu, pN = pN, pK = pK, nExp = nExp_poi, reuse.pANN = FALSE)

  ## store results
  colData(sce)$dblFinder_score = seu@meta.data[, paste0("pANN_", cn)]
  colData(sce)$dblFinder_call  = seu@meta.data[, paste0("DF.classifications_", cn)]
  metadata(sce)$doubletfinder = list("params" = cn,
                                     "bcmvn" = bcmvn,
                                     "pK" = pK)

  ## run with prior cluster annotation information
  ## i. e., if user-supplied literature supported cell type annotations available
  if( !is.null(sce$ANA_cluster_labels) & use_celltype_anno){
    ## Homotypic Doublet Proportion Estimate
    annotations <- sce$ANA_cluster_labels
    homotypic.prop <- modelHomotypic(annotations)
    homotypic.prop <- modelHomotypic(annotations)
    nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
    cn.adj <- paste(pN, pK, nExp_poi.adj, sep = "_")
    cat(cn.adj, "\n")

    ## "Generate new list of doublet classificatons from existing pANN vector to save time"
    ## cn has pANN values calc using Pois assumption
    seu <- doubletFinder(seu, pN = pN, pK = pK, nExp = nExp_poi.adj, reuse.pANN = paste0("pANN_", cn))
    ## with higher exp no. of doublets
    df_hi.lo <- seu@meta.data[, paste0("DF.classifications_", cn)]
    ## with adj for homotypic, so lower no. of expected doublets
    df.adj <- seu@meta.data[, paste0("DF.classifications_", cn.adj)]
    ## set those that show up in adj one as high confidence, rest to low
    df_hi.lo[df_hi.lo == "Doublet" & df.adj == "Singlet"] <- "Doublet_lo"
    df_hi.lo[df_hi.lo == "Doublet"] <- "Doublet_hi"

    ## store results
    colData(sce)$doubletfinder_call_hilo  = df_hi.lo
    metadata(sce)$doubletfinder_call_hilo = list("params" = cn.adj)
  }
  return(sce)
}
