

mkLB <- function(resp,pred){
  res = rep("TN",length(pred))
  res[resp & pred] = "TP"
  res[resp & (!pred)] = "FN"
  res[(!resp) & pred] = "FP"
  return(res)
}

mkTile <- function(nme,dat){
  mclrs = c(rgb(.8,0,0,.5),rgb(0,0,0,.7),rgb(0,0,.8,.5),rgb(0,0,0,.7))
  names(mclrs) = c("FP","TP","FN","TN")
  mcs  <- scale_fill_manual(name = nme, values = mclrs)
  p = ggplot(dat,aes_string(x=nme,y="libsize",fill=nme)) + geom_violin(col=NA,scale="area") + theme(legend.position="none") + mcs
  return(p)
}

plotLibSzeStrat <- function(sce,filename){

  dat = data.frame(call = sce$dbl_anno)

  dat$annotation  = as.integer(dat$call)
  frc               = mean(dat$call)
  dat$cxds          = binFu(sce$cxds_score,frc)
  dat$scrublet      = binFu(sce$scrublet_score,frc)
  dat$bcds          = binFu(sce$bcds_score,frc)
  dat$dblFinder     = binFu(sce$dblFinder_score,frc)
  dat$dblDetection  = binFu(sce$dblDetection_score,frc)
  dat$hybrid        = binFu(sce$hybrid_score,frc)
  dat$features      = binFu(sce$numfeat_score,frc)
  dat$dblCells      = binFu(sce$dblCells_score,frc)
  dat$cxds          = mkLB(dat$call,dat$cxds)
  dat$bcds          = mkLB(dat$call,dat$bcds)
  dat$scrublet      = mkLB(dat$call,dat$scrublet)
  dat$dblFinder     = mkLB(dat$call,dat$dblFinder)
  dat$dblCells      = mkLB(dat$call,dat$dblCells)
  dat$dblDetection  = mkLB(dat$call,dat$dblDetection)
  dat$hybrid        = mkLB(dat$call,dat$hybrid)
  dat$libsize       = log10(Matrix::colSums(counts(sce)))

  pdf("/dev/null")
  plst = lapply(c("cxds","bcds","hybrid","dblDetection","dblFinder","dblCells","scrublet"),
                mkTile,dat)
  dev.off()
  pdf(file=filename,width=14,height=2)
  print(plot_grid(plotlist=plst,nrow=1))
  dev.off()
}
