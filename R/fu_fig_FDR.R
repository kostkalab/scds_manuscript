
library(Matrix)
library(SingleCellExperiment)
library(knitr)
library(kableExtra)
library(scds)

sce.chcl = readRDS("docker_results_scds_fin_perf/results/sce_chcl.rds")
sce.chpb = readRDS("docker_results_scds_fin_perf/results/sce_chpb.rds")
sce.demu = readRDS("docker_results_scds_fin_perf/results/sce_demu.rds")
sce.hgmm = readRDS("docker_results_scds_fin_perf/results/sce_hgmm.rds")

#- get the number of false false positives for each score
get_nfps <- function(scr,ann,method){
  ord   = order(scr,decreasing=TRUE)
  scr_s = scr[ord]
  ann_s = ann[ord]
  nfp_s = cumsum(!ann_s)
  ntp_s = cumsum( ann_s)
  fdr_s = nfp_s/(nfp_s+ntp_s)

  res   = cbind(seq_len(length(ord)),scr_s,ann_s,nfp_s,ntp_s,fdr_s)
  res   = data.frame(res)
  res$method = method
  res$ftp    = res$ntp/max(res$ntp)
  colnames(res) = c("rank","scr","ann","nfp","ntp","fdr","method","ftp")
  rownames(res) = NULL
  return(res)
}

afu <- function(meth,sce,short=TRUE){
  scr = colData(sce)[,paste(meth,"score",sep="_")]
  ann = sce$dbl_anno

  res = get_nfps(scr,ann,meth)
  if(short) res = res[seq_len(sum(ann)),]
  return(res)
}

#- data sets and parameter ranges
dats <- c("chcl", "chpb", "demu", "hgmm")
mets <- c("dblCells","dblFinder","dblDetection","scrublet","cxds","bcds","hybrid", "libsize","numfeat")

dfu <- function(sce,dname,short=TRUE){
  tmp  = sapply(mets,afu,sce,short=short, simplify=FALSE)
  df   = tmp[[1]]
  for(i in 2:length(tmp)) df = rbind(df, tmp[[i]])
  df$data = dname
  return(df)
}
