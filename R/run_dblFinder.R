library(scran)
library(Seurat)
library(DoubletFinder)


## UPDATE: to skip for NaN values, which occur in micm data analysis
## use this instead of DoubletFinder::summarizeSweep
summarizeSweep_NA <- function (sweep.list, GT = FALSE, GT.calls = NULL) {
    require(KernSmooth)
    require(ROCR)
    require(modes)
    pK <- c(5e-04, 0.001, 0.005, seq(0.01, 0.3, by = 0.01))
    pN <- seq(0.05, 0.3, by = 0.05)
    if (GT == TRUE) {
        sweep.stats <- as.data.frame(matrix(0L, nrow = length(sweep.list),
                                            ncol = 4))
        colnames(sweep.stats) <- c("pN", "pK", "AUC", "BCreal")
        sweep.stats$pN <- factor(rep(pN, each = length(pK), levels = pN))
        sweep.stats$pK <- factor(rep(pK, length(pN), levels = pK))
    }
    if (GT == FALSE) {
        sweep.stats <- as.data.frame(matrix(0L, nrow = length(sweep.list),
                                            ncol = 3))
        colnames(sweep.stats) <- c("pN", "pK", "BCreal")
        sweep.stats$pN <- factor(rep(pN, each = length(pK), levels = pN))
        sweep.stats$pK <- factor(rep(pK, length(pN), levels = pK))
    }
    for (i in 1:length(sweep.list)) {
        res.temp <- sweep.list[[i]]
        ## Update: handle NaN cases: SKIP
        if(any(is.na(res.temp))){
            message("Skipped: ", names(sweep.list)[i], " (it is NA).")
            next
        }
        gkde <- approxfun(bkde(res.temp$pANN, kernel = "normal"))
        x <- seq(from = min(res.temp$pANN), to = max(res.temp$pANN),
            length.out = nrow(res.temp))
        sweep.stats$BCreal[i] <- bimodality_coefficient(gkde(x))
        if (GT == FALSE) {
            next
        }
        meta <- as.data.frame(matrix(0L, nrow = nrow(res.temp),
            ncol = 2))
        meta[, 1] <- GT.calls
        meta[, 2] <- res.temp$pANN
        train.ind <- sample(1:nrow(meta), round(nrow(meta)/2),
            replace = FALSE)
        test.ind <- (1:nrow(meta))[-train.ind]
        colnames(meta) <- c("SinDub", "pANN")
        meta$SinDub <- factor(meta$SinDub, levels = c("Doublet",
            "Singlet"))
        model.lm <- glm(SinDub ~ pANN, family = binomial(link = "logit"),
            data = meta, subset = train.ind)
        prob <- predict(model.lm, newdata = meta[test.ind, ],
            type = "response")
        ROCpred <- ROCR::prediction(predictions = prob, labels = meta$SinDub[test.ind])
        perf.auc <- ROCR::performance(ROCpred, measure = "auc")
        sweep.stats$AUC[i] <- perf.auc@y.values[[1]]
    }
    return(sweep.stats)
}


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
  sweep.stats <- summarizeSweep_NA(sweep.res.list) ## GT = FALSE is default

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
  colData(sce)$dblFinder_call  = seu@meta.data[, paste0("DF.classifications_", cn)] == "Doublet"
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
