
library(pROC)
library(PRROC)

scorePred <- function(resp, pred){
#=================================

  #remove NAs if present; probably shoud tally
  resp = resp[!is.na(pred)]
  pred = pred[!is.na(pred)]

  pr  = pr.curve(scores.class0 = pred, weights.class0 = resp)

  tmp = coords(roc(response = resp, predictor=pred),x=quantile(pred,.95), input="threshold",ret=c("precision","recall"))
  a.0 = roc(response=resp, predictor=pred)$auc
  a.1 = roc(response=resp, predictor=pred,partial.auc=c(1,.9)  ,partial.auc.focus="spe",partial.auc.correct=TRUE)$auc
  a.2 = roc(response=resp, predictor=pred,partial.auc=c(1,.95) ,partial.auc.focus="spe",partial.auc.correct=TRUE)$auc
  a.3 = roc(response=resp, predictor=pred,partial.auc=c(1,.975),partial.auc.focus="spe",partial.auc.correct=TRUE)$auc
  a.4 = roc(response=resp, predictor=pred,partial.auc=c(1,.99) ,partial.auc.focus="spe",partial.auc.correct=TRUE)$auc
  p.1 = tmp[1]
  r.1 = tmp[2]
  ap  = pr$auc.davis.goadrich

  rn  = c("AUC","pAUC900","pAUC950","pAUC975","pAUC990","prec05","rec05","AUPRC")
  res = c(a.0, a.1, a.2, a.3, a.4, p.1, r.1, ap)
  names(res) = rn
  return(res)
}

mk_perf_tab <-function(sce,ord=TRUE){
#====================================

  r1  = scorePred(resp=sce$dbl_anno,pred=sce$cxds_score)
  r2  = scorePred(resp=sce$dbl_anno,pred=sce$bcds_score)
  r3  = scorePred(resp=sce$dbl_anno,pred=sce$hybrid_score)
  r4  = scorePred(resp=sce$dbl_anno,pred=sce$libsize_score)
  r5  = scorePred(resp=sce$dbl_anno,pred=sce$numfeat_score)
  r6  = scorePred(resp=sce$dbl_anno,pred=sce$scrublet_score)
  r7  = scorePred(resp=sce$dbl_anno,pred=sce$dblDetection_score)
  r8  = scorePred(resp=sce$dbl_anno,pred=sce$dblFinder_score)
  r9  = scorePred(resp=sce$dbl_anno,pred=sce$dblCells_score)

  tab           = rbind(r1,r2,r3,r4,r5,r6,r7,r8,r9)
  rownames(tab) = c("cxds","bcds","hybrid","libsize","features",
                    "scrublet","dblDetection","dblFinder","dblCells")
  tab           = round(tab,2)
  oo            = 1:nrow(tab)
  if(ord) oo    = order(tab[,1],tab[,8])
  return(tab[oo,])
}
