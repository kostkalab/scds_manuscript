
library(sctransform)
library(xgboost)

bcds_nn <- function(sce, ntop=500, srat=1, verb=FALSE, retRes=FALSE,
                 nmax="tune", varImp=FALSE, estNdbl = FALSE){

  #- check first argument (sce)
  if(!is(sce,"SingleCellExperiment"))
      stop('First argument (sce) needs to be of class SingleCellExperiment')
  if( !("counts" %in% names(assays(sce))) )
      stop('First argument (sce) needs to have a "counts" assay')
  if(!is(counts(sce),"sparseMatrix"))
      counts(sce) = Matrix::Matrix(counts(sce),sparse=TRUE)

  #- select variable genes
  #-----------------------
  if(verb) message("-> normalizing and selecting genes\n")
  sce_data <- assay(sce, "counts")
  future::plan(strategy = 'multicore', workers = 4)
  options(future.globals.maxSize = 10 * 1024 ^ 3)
  set.seed(100)
  vst_out <- sctransform::vst(sce_data,
                                latent_var = c('log_umi'),
                                return_gene_attr = TRUE,
                                return_cell_attr = TRUE,
                                show_progress = TRUE,
                                min_cells = 0.01 * ncol(sce_data))
  vrs <- vst_out$gene_attr$residual_variance
  names(vrs) <- rownames(vst_out$gene_attr)
  hvg  <- order(vrs, decreasing = TRUE)[seq_len(ntop)]
  hvg.genes <- rownames(vst_out$gene_attr)[hvg]

  #- Add artificial doublets (joint mormalization)
  #-------------------------------------------
  ind1 = Matrix::rowSums(counts(sce) > 0) > 0.01 * ncol(sce)
  if(verb) message("-> simulating doublets\n")
  p1  = sample(seq_len(ncol(sce)), srat*ncol(sce), replace = TRUE)
  p2  = sample(seq_len(ncol(sce)), srat*ncol(sce), replace = TRUE)
  ## simulated data
  sim_data <- counts(sce)[, p1] + counts(sce)[, p2]
  colnames(sim_data) <- paste0("sim_", seq_len(ncol(sim_data)))

  ## artificial dataset with simulated doublets included for the original cells
  new_data <- cbind(sce_data, sim_data)[ind1, ]
  ## normalize
  set.seed(100)
  vst_out1 <- sctransform::vst(new_data,
                              latent_var = c('log_umi'),
                              return_gene_attr = TRUE,
                              return_cell_att = TRUE,
                              show_progress = TRUE)
  norm_data <- vst_out1$y
  X = t(norm_data[hvg.genes[hvg.genes %in% rownames(norm_data)], ]) #- some m ight not be there`
  #- learn classifier
  if(verb) message("-> training classifier\n")
  colnames(X) = paste("GEN",seq_len(ncol(X)),sep="_") #- ranger picky
  y           = c(rep(-1,ncol(sce)),rep(1,length(p1)))

  #======================================
  #- NOTE: rest is same as "normal" bcds
  #======================================
  mm  = xgb.DMatrix(X,label=(y+1)/2)
  #- fixed rounds:
  if(nmax != "tune"){
    res    = NA
    varImp = FALSE
    retRes = FALSE
    pre    = xgboost(mm,nrounds=nmax,tree_method="hist",
                     nthread = 2, early_stopping_rounds = 2, subsample=0.5,
                     objective = "binary:logistic",verbose=0)
    sce$bcds_score = stats::predict(pre, newdat= mm[seq_len(ncol(sce)),])
    #- get doublet calls
    if(estNdbl){
      dbl.pre = stats::predict(pre, newdat= mm[seq(ncol(sce)+1,nrow(X)),])
      est_dbl = get_dblCalls_ALL(sce$bcds_score,dbl.pre,rel_loss=srat)
      if(is.null(metadata(sce)$cxds)) metadata(sce)$cxds = list()
      metadata(sce)$bcds$ndbl = est_dbl
      metadata(sce)$bcds$sim_scores = dbl.pre
      sce$bcds_call = sce$bcds_score >= est_dbl["balanced","threshold"]
    }
  #- learning rounds with CV:
  } else {
    res = xgb.cv(data =mm, nthread = 2, nrounds = 500, objective = "binary:logistic",
                 nfold=5,metrics=list("error"),prediction=TRUE,
                 early_stopping_rounds=2, tree_method="hist",subsample=0.5,verbose=0)
    ni  = res$best_iteration
    ac  = res$evaluation_log$test_error_mean[ni] + 1*res$evaluation_log$test_error_std[ni]
    ni  = min(which( res$evaluation_log$test_error_mean <= ac  ))
    nmax = ni
    pre = xgboost(mm,nrounds=nmax,tree_method="hist",
                  nthread = 2, early_stopping_rounds = 2, subsample=0.5,
                  objective = "binary:logistic",verbose=0)
    sce$bcds_score = res$pred[seq_len(ncol(sce))]
    #- get doublet calls
    if(estNdbl){
      dbl.pre = stats::predict(pre, newdat= mm[seq(ncol(sce)+1,nrow(X)),])
      est_dbl = get_dblCalls_ALL(sce$bcds_score,dbl.pre,rel_loss=srat)
      if(is.null(metadata(sce)$cxds)) metadata(sce)$cxds = list()
      metadata(sce)$bcds$ndbl = est_dbl
      metadata(sce)$bcds$sim_scores = dbl.pre
      sce$bcds_call = sce$bcds_score >= est_dbl["balanced","threshold"]

    }
  }
  #- variable importance
  if(varImp){
    if(verb) message("-> calculating variable importance\n")
    vimp = xgb.importance(model=pre)
    vimp$col_index = match(vimp$Feature,colnames(X))
  }
  #- result
  if(retRes){
    hvg_bool                    = (seq_len(nrow(sce))) %in% which(ind1)
    hvg_bool[hvg_bool]          = (seq_len(sum(hvg_bool))) %in% hvg
    hvg_ord                     = rep(NA,nrow(sce))
    hvg_ord[hvg_bool]           = which(ind1)[hvg]
    rowData(sce)$bcds_hvg_bool  = hvg_bool
    rowData(sce)$bcds_hvg_ordr  = hvg_ord
    metadata(sce)$bcds_res_cv   = res
    metadata(sce)$bcds_res_all  = pre
    metadata(sce)$bcds_nmax     = nmax
    if(varImp){
      vimp$gene_index             = hvg_ord[hvg_bool][vimp$col_index]
      metadata(sce)$bcds_vimp     = vimp[seq_len(100),-c(1,5)]
    }
  }

  if(verb) message("-> done.\n\n")
  return(sce)

}
