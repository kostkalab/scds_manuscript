library(pamr)
library(scds)
library(SingleCellExperiment)
library(xgboost)
library(pROC)

getDAT <- function(sce){
#=======================

  snms = c(  "cxds_score",
             "bcds_score",
             "libsize_score",
             "libsize_score",
             "dblCells_score" ,
             "dblDetection_score" ,
             "dblFinder_score",
             "scrublet_score","dblCells_score")

  X    = colData(sce)[,snms]
  ii   = which(rowSums(is.na(X))>0)
  if(length(ii)>0) {
    X    = X[-ii,]
    y    = as.factor(as.integer(sce$dbl_anno)[-ii])
  } else {
    y    = as.factor(as.integer(sce$dbl_anno))
  }
  #- use ranks so we can learn -> predict on different sces
  X    = apply(X,2,function(x) {x = rank(x)/length(x); return(x)})
  X    = data.frame(X)

  return(list(X=X,y=y))
}

getModel <- function(X,y){
#=========================
  mm   = xgb.DMatrix(as(X,"matrix"),label=as.numeric(y)-1)
  res  = xgb.cv(data =mm, nthread = 2, nrounds = 100, objective = "binary:logistic",
                nfold=5,metrics=list("auc"),prediction=FALSE,
                early_stopping_rounds=7, scale_pos_weight=sum(y==0)/sum(y==1))
  ni   = res$best_iteration
  nmax = ni
  pre  = xgboost(mm,nrounds=nmax,metrics=list("auc"),
                 nthread = 2, early_stopping_rounds = 7,
                 objective = "binary:logistic",)
  return(pre)
}

assessModel <- function(mod,X,y){
#================================

  mat = xgb.DMatrix(as(X,"matrix"))
  pr  = predict(mod,mat)
  rc  = roc(response = y, predictor = pr)
  return(rc$auc)

}

assessModelCV <- function(X,y,nf=5){
#===================================

  bf = pamr:::balanced.folds(y,nf)

  pr = NULL
  yr = NULL
  for(i in seq_len(nf)){
    X.train = X[-sort(bf[[i]]),]
    y.train = y[-sort(bf[[i]])]
    X.test  = X[bf[[i]],]
    yr      = c(yr,y[bf[[i]]])
    md      = getModel(X.train,y.train)
    mat     = xgb.DMatrix(as(X.test,"matrix"))
    pr      = c(pr, predict(md,mat))
  }

  rc = roc(response=yr,predictor = pr)
  return(rc$auc)

}

assessBestFeature <- function(dat){
#==================================

  aucs = sapply(colnames(dat$X), function(cn) roc(resp=dat$y,pred=dat$X[,cn])$auc)
  return(max(aucs))
}
