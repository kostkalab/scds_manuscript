library(scran)
library(pROC)
library(Matrix)
library(reticulate)
use_python("/opt/user/rstudio/miniconda/bin/python", required = FALSE)
use_condaenv("base",conda="/opt/user/rstudio/miniconda/bin/conda")

dbd        = import("doubletdetection")

run_dblDetection <- function(sce){
  ddclf      = dbd$BoostClassifier()
  res.dbd    = ddclf$fit(t(counts(sce)))
  preds      = res.dbd$predict()
  scores     = colMeans(res.dbd$all_log_p_values_, na.rm = TRUE)

  colData(sce)$dblDetection_score = -scores
  colData(sce)$dblDetection_call  = preds
  return(sce)
}

