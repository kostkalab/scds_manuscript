library(scran)
library(scater)

run_dblCells <- function(sce){
  cat("normalize ... \n")
  set.seed(1000)
  clusters <- quickCluster(sce, method = "igraph", min.mean = 0.1)
  table(clusters)
  sce <- computeSumFactors(sce, clusters = clusters, min.mean = 0.1)
  summary(sizeFactors(sce))
  sce <- normalize(sce)
  assayNames(sce)

  cat("remove technical noise ... \n")
  tech.trend <- makeTechTrend(x = sce)
  fit <- trendVar(sce, use.spikes = FALSE)

  set.seed(12345)
  sce <- denoisePCA(sce, technical = tech.trend, approximate = TRUE)
  ncol(reducedDim(sce))

  cat("run doubletCells... \n")
  set.seed(100)
  dbl.dens <- doubletCells(sce, approximate=TRUE)
  summary(dbl.dens)

  colData(sce)$dblCells_score = dbl.dens
  colData(sce)$dblCells_call  = rep(NA,ncol(sce))
  return(sce)
}

