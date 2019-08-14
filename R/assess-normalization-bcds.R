
library(SingleCellExperiment)
library(scds)
library(pROC)
library(PRROC)
library(knitr)
library(kableExtra)

sce.chcl = readRDS("docker_results_scds_fin_perf/results/sce_chcl.rds")
sce.chpb = readRDS("docker_results_scds_fin_perf/results/sce_chpb.rds")
sce.demu = readRDS("docker_results_scds_fin_perf/results/sce_demu.rds")
sce.hgmm = readRDS("docker_results_scds_fin_perf/results/sce_hgmm.rds")

source("./R/fu_tab_norm_bcds.R")
source("./R/fu_tab_performance.R")

comp <- function(sce, nme){
  s1 = bcds(sce,verb=T)
  s2 = bcds_nn(sce,verb=T)
  p1 = scorePred(s1$dbl_anno,s1$bcds_score)
  p2 = scorePred(s2$dbl_anno,s2$bcds_score)
  res = cbind(p1,p2)
  colnames(res) = c(paste(nme,"bcds",sep="_"),paste(nme,"bcds-new",sep="_"))
  return(round(res,3))
}

t1 = comp(sce.chcl,"chcl")
t2 = comp(sce.chpb,"chpb")
t3 = comp(sce.demu,"demu")
t4 = comp(sce.hgmm,"hgmm")

tab = cbind(t1,t2,t3,t4)

tabres <- t(tab) %>% kable(format="latex",booktabs=TRUE) %>%
kable_styling()                                          %>%
group_rows("chcl",1,2)                                   %>%
group_rows("chpb",3,4)                                   %>%
group_rows("demu",5,6)                                   %>%
group_rows("hgmm",7,8)


#- NEW a bit better for demu , chpb but not the others.  small improvement, not worth the additonal running times. Convergence issues with the new norm.
