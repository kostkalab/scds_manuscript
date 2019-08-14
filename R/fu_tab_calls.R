

#- precision
get_pre = function(ann,pred){
#----------------------------

  TP = sum(ann  & pred)
  FP = sum(pred & (!ann))
  return(TP/(TP+FP))
}


#- number of false neagives (doublets missed)
get_FN = function(ann,pred){
#----------------------------
  sum(ann & (!pred))
}

#- recall
get_rec = function(ann,pred){
#----------------------------

  TP = sum(ann & pred)
  FN = sum(ann & (!pred))
  return(TP/(TP+FN))
}

#- balanced accuracy
get_bac = function(ann,pred){
#----------------------------

  TP  = sum(ann    & pred)
  FN  = sum(ann    & (!pred))
  TN  = sum((!ann) & (!pred))
  FP  = sum((!ann) & pred)
  TPR = TP/(TP+FN)
  TNR = TN/(TN+FP)
  return((TPR+TNR)/2)
}

mk_call_tab <- function(sce){
#----------------------------

  rowfu = function(x){
    ind = !is.na(x)
    x   = x[ind]
    y   = sce$dbl_anno[ind]
    c(sum(x),get_pre(y,x),get_rec(y,x),get_bac(y,x), get_FN(y,x))
  }

  r1   = rowfu(sce$dblFinder_call)
  r2   = rowfu(sce$dblDetection_call)
  r3   = rowfu(sce$scrublet_call)
  r4   = rowfu(sce$dblDecon_call)
  r5   = rowfu(sce$cxds_call) #- cxds - balanced
#  r6   = rowfu(sce$cxds_score>=quantile(metadata(sce)$cxds$sim_scores,0.01))
#  r7   = rowfu(sce$cxds_score>=quantile(metadata(sce)$cxds$sim_scores,0.1))
  r8   = rowfu(sce$bcds_call) #- cxds - balanced
#  r9   = rowfu(sce$bcds_score>=quantile(metadata(sce)$bcds$sim_scores,0.01))
#  r10  = rowfu(sce$bcds_score>=quantile(metadata(sce)$bcds$sim_scores,0.1))
  r11  = rowfu(sce$hybrid_call) #- cxds - balanced
#  r12  = rowfu(sce$hybrid_score>=quantile(metadata(sce)$hybrid$sim_scores,0.01))
#  r13  = rowfu(sce$hybrid_score>=quantile(metadata(sce)$hybrid$sim_scores,0.1))

#  res = rbind(r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,r13)
  res = rbind(r1,r2,r3,r4,r5,r8,r11)
  colnames(res) = c("num_dbl","precision","recall","balAcc","FN")
  res = data.frame(round(res,2))

  #res$cut_type=(c(rep("built_in",4),rep(c("balanced","0.01","0.1"),3)))
  res$cut_type=(c(rep("built_in",4),rep(c("balanced"),3)))

#  res$method = (c("dblFinder","dblDetection","scrublet","dblDecon",rep("cxds",3),rep("bcds",3),rep("hybrid",3)))
  res$method = c("dblFinder", "dblDetection","scrublet","dblDecon","cxds","bcds","hybrid")


  return(res)
}
