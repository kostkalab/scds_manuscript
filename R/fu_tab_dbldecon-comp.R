
countTPFP <- function(resp,pred){
#================================

  res        = rep(NA,4)
  names(res) = c("TP","FP","FN","TN")
  res["TP"]  = sum(resp & pred)
  res["FP"]  = sum((!resp) & pred)
  res["TN"]  = sum((!resp) & (!pred))
  res["FN"]  = sum(resp & (!pred))
  return(res)
}

binFu <- function(x,frc){
#========================

  numpos = round(length(x)*frc)
  xin    = order(x,decreasing=TRUE)[1:numpos]
  res    = as.integer(1:length(x) %in% xin)
  return(res)
}

mk_dd_tab <- function(sce){
  #==========================
  ## maybe change dblDecon calls to logical in orig func
  if(!is.logical(sce$dblDecon_call)){
      dcall = sce$dblDecon_call
      do <- rep(FALSE, ncol(sce))
      do[which(dcall == "Doublet")] = TRUE
      names(do) <- colnames(sce)
      sce$dblDecon_call <- do
  }
  frac = mean(sce$dblDecon_call)
  res = rbind(
    countTPFP(sce$dbl_anno,sce$dblDecon_call)                 ,
    countTPFP(sce$dbl_anno,binFu(sce$libsize_score,frac))     ,
    countTPFP(sce$dbl_anno,binFu(sce$numfeat_score,frac))     ,
    countTPFP(sce$dbl_anno,binFu(sce$dblCells_score,frac))    ,
    countTPFP(sce$dbl_anno,binFu(sce$scrublet_score,frac))    ,
    countTPFP(sce$dbl_anno,binFu(sce$dblFinder_score,frac))   ,
    countTPFP(sce$dbl_anno,binFu(sce$dblDetection_score,frac)),
    countTPFP(sce$dbl_anno,binFu(sce$cxds_score,frac))        ,
    countTPFP(sce$dbl_anno,binFu(sce$bcds_score,frac))        ,
    countTPFP(sce$dbl_anno,binFu(sce$hybrid_score,frac))
  )
  rownames(res) = c("dblDecon","libsize","features","dblCells","scrublet",
                    "dblFinder","dblDetection","cxds","bcds","hybrid")
  res = as.data.frame(res)
  res$sen = res[,"TP"]/(res[,"TP"]+res[,"FN"])
  res$spe = res[,"TN"]/(res[,"TN"]+res[,"FP"])
  res$pre = res[,"TP"]/(res[,"TP"]+res[,"FP"])
  return(res)
}
