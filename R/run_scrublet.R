library(scran)
library(pROC)
library(Matrix)
library(reticulate)
library(reticulate)
use_python("/opt/user/rstudio/miniconda/bin/python", required = FALSE)
use_condaenv("base",conda="/opt/user/rstudio/miniconda/bin/conda")

scr        = import("scrublet")

run_scrublet <- function(sce){
#=============================
  scrub      = scr$scrublet$Scrublet(t(counts(sce)))
  res        = scrub$scrub_doublets()
  names(res) = c("dbl_score","dbl_call")
  colData(sce)$scrublet_score = res$dbl_score
  colData(sce)$scrublet_call  = res$dbl_call
  return(sce)
}

