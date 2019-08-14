#- START TBR
library(ggplot2)
library(dplyr)
library(kableExtra)
library(SingleCellExperiment)
library(pROC)
library(PRROC)
library(ggbeeswarm)

sce.chcl = readRDS("./results/sce_chcl.rds")
sce.chpb = readRDS("./results/sce_chpb.rds")
sce.demu = readRDS("./results/sce_demu.rds")
sce.hgmm = readRDS("./results/sce_hgmm.rds")

#- need performance measures
#- performance; see ./R/fu_tab_performance.R
get_auroc <- function(resp,pred) roc(response=resp, predictor=pred)$auc
get_auprc <- function(resp,pred) pr.curve(scores.class0 = pred, weights.class0 = resp)$auc.davis.goadrich

set.seed(23124)
NCEL = 2000
NREP = 100

#- data sets and parameter ranges
dats <- c("chcl", "chpb", "demu", "hgmm")
mets <- c("dblCells","dblFinder","dblDetection","scrublet","cxds","bcds","hybrid")

method_names   = paste("run_",mets,sep="")
method_scripts = paste("./R/run_",mets,".R",sep="")
for(i in seq_along(method_scripts)) source(method_scripts[i])

RES = NULL
for(dat in dats){
  message(dat)
  #- read in data set
  data_file     = paste("./data/",dat,"/proc/sce_",dat,".rds",sep="")
  sce           = readRDS(data_file)

  afu <- function(sce){
  #--------------------

    sce.rs = sce[,sample.int(ncol(sce),size=NCEL)]
    sce.rs = sce.rs[Matrix::rowSums(counts(sce.rs) > 0) > 0,]
    res           = matrix(NA,ncol=2,nrow=length(method_names)+2)
    rownames(res) = c(mets,"libsize","features")
    colnames(res) = c("AUROC","AUPRC")

    #- run methods on sub-sampled data
    for(i in seq_along(method_names)){
      sce.rs         = do.call(method_names[i],list(sce=sce.rs))
      met            = sub("run_","",method_names[i])
      scr            = colData(sce.rs)[,paste(met,"score",sep="_")]
      ind            = !is.na(scr) #- rarely happens, but still
      res[i,"AUROC"] = get_auroc(resp=sce.rs$dbl_anno[ind], pred=scr[ind])
      res[i,"AUPRC"] = get_auprc(resp=sce.rs$dbl_anno[ind], pred=scr[ind])
    }

    #- baseline methods
    scr                       = Matrix::colSums(counts(sce.rs))
    res[nrow(res)-1,"AUROC"]  = get_auroc(resp=sce.rs$dbl_anno, pred=scr)
    res[nrow(res)-1,"AUPRC"]  = get_auprc(resp=sce.rs$dbl_anno, pred=scr)
    scr                       = Matrix::colSums(counts(sce.rs)>0)
    res[nrow(res),"AUROC"]    = get_auroc(resp=sce.rs$dbl_anno, pred=scr)
    res[nrow(res),"AUPRC"]    = get_auprc(resp=sce.rs$dbl_anno, pred=scr)
    return(res)
  }
  res = replicate(NREP,afu(sce))
  res = lapply(seq_len(dim(res)[3]), function(i) { tmp=as.data.frame(res[,,i]); tmp$rs_id=i; tmp })
  res = lapply(res,function(x) {x$data = dat; x})
  print(res)
  RES = c(RES,res)
  rm(res)
}

#- flatten list to data frame
#-----------------------------
pRES = lapply(RES,function(x) {x$method=rownames(x); rownames(x) = NULL; x})
df   = pRES[[1]]
for(i in seq_along(2:length(pRES))) df = rbind(df, pRES[[i]])

#- save temp results
saveRDS(df,"./results/tmp-resamp.rds")

#- PLOTS
#=======
theme_set(theme_grey())
mk_plot <- function(df, ystr){
#-----------------------------
  p <- ggplot(df, aes_string(x="method",y=ystr,col="method")) +
       geom_quasirandom(cex= 0.2) +
       stat_summary(fun.data = mean_sdl, color = "black", size = 0.2) +
       facet_wrap(~data) +
       xlab("method") + ylab(ystr) +
       theme(axis.text.x  = element_text(size = 10, angle = 90,
          vjust=0.5), axis.text.y = element_text(size = 12),
          strip.text.x = element_text(size = 12),
          strip.background = element_blank(), legend.position = "none") +
       guides(colour = guide_legend(override.aes = list(size = 2),
          title = "method"))
  return(p)
}

df$method <- df$method %>% as.factor %>% reorder(df$AUROC)
p1 <- mk_plot(df,"AUROC")

#- plot AUPRC
df$method <- df$method %>% as.factor %>% reorder(df$AUPRC)
p2 <- mk_plot(df,"AUPRC")

#- save to files:
jpeg(filename = "./results/fig_resampling_auroc.jpeg", width = 6, height = 7, units = "in", res = 1200)
print(p1)
dev.off()
jpeg(filename = "./results/fig_resampling_auprc.jpeg", width = 6, height = 7, units = "in", res = 1200)
print(p2)
dev.off()
