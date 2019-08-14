
library(SingleCellExperiment)
library(scds)
library(scater)
library(cowplot)
library(ggplot2)
library(reshape2)
theme_set(theme_grey())
library(dplyr)

sce.chcl <- readRDS("./results/sce_chcl.rds") %>% calculateQCMetrics
sce.chpb <- readRDS("./results/sce_chpb.rds") %>% calculateQCMetrics
sce.demu <- readRDS("./results/sce_demu.rds") %>% calculateQCMetrics
sce.hgmm <- readRDS("./results/sce_hgmm.rds") %>% calculateQCMetrics

mk_plot <- function(sce,nme=""){
#===============================

    #- make the c
    ndbl  = sum(sce$dbl_anno)

    gc    = function(x,num) {
      xin    = order(x,decreasing=TRUE)[1:num]
      res    = as.integer(1:length(x) %in% xin)
      return(res)
    }
    calls = data.frame(
      bcds_call         = gc(sce$bcds_score,ndbl),
      cxds_call         = gc(sce$cxds_score,ndbl),
      hybrid_call       = gc(sce$hybrid_score,ndbl),
      dblDetection_call = gc(sce$dblDetection_score,ndbl),
      dblFinder_call    = gc(sce$dblFinder_score,ndbl),
      dblCells_call     = gc(sce$dblCells_score,ndbl),
      libsize_call      = gc(sce$libsize_score,ndbl),
      numfeat_call      = gc(sce$numfeat_score,ndbl),
      scrublet_call     = gc(sce$scrublet_score,ndbl),
      anno              = sce$dbl_anno
    )


    df = colData(sce)[,c("log10_total_counts",
              "pct_counts_in_top_50_features")]
    df$all_missed = (rowSums(calls[,1:9]) == 0) & calls[,10]
    df = df[sce$dbl_anno,]
    df = df[rowSums(is.na(df))==0,]
    df = as.data.frame(df)
    xx = melt(df)
    p = ggplot(xx,aes(x=all_missed,y=value)) + geom_violin(fill="black") + facet_wrap(~variable,scales="free") + ylab("") + ggtitle(nme)
    return(p)
}

#- doublets missed by all methods
#================================

p.chcl = mk_plot(sce.chcl,nme="chcl")
p.chpb = mk_plot(sce.chpb,nme="chpb")
p.demu = mk_plot(sce.demu,nme="demu")

pdf(file="./results/fig_dbl-EDA.pdf",width=3*3*1.5,height=1*3*1.5)
plot_grid(p.chcl,p.chpb,p.demu,nrow=1)
dev.off()



#- scrublet on chcl
#==================

sce = sce.chcl
ndbl  = sum(sce$dbl_anno)

gc    = function(x,num) {
  xin    = order(x,decreasing=TRUE)[1:num]
  res    = as.integer(1:length(x) %in% xin)
  return(res)
}

calls = data.frame(
  bcds_call         = gc(sce$bcds_score,ndbl),
  cxds_call         = gc(sce$cxds_score,ndbl),
  hybrid_call       = gc(sce$hybrid_score,ndbl),
  dblDetection_call = gc(sce$dblDetection_score,ndbl),
  dblFinder_call    = gc(sce$dblFinder_score,ndbl),
  dblCells_call     = gc(sce$dblCells_score,ndbl),
  libsize_call      = gc(sce$libsize_score,ndbl),
  numfeat_call      = gc(sce$numfeat_score,ndbl),
  scrublet_call     = gc(sce$scrublet_score,ndbl),
  anno              = as.integer(sce$dbl_anno)
)

scrublet_only = apply(calls,1,function(x) all(x == c(0,0,0,0,0,0,0,0,1,1)))

sum(scrublet_only) #- 29
anymeth = apply(calls,1,function(x) x[10]==1 & any(x[-10]==1) )
anymeth[scrublet_only] = FALSE
sum(sce$log10_total_counts[scrublet_only] < min(sce$log10_total_counts[anymeth])) #- 6


calls_tp = calls * calls[,7]    ; calls_tp = calls_tp[,-7]
calls_fp = calls * (!calls[,7]) ; calls_fp = calls_fp[,-7]
colnames(calls_tp) = paste(colnames(calls_tp),"tp",sep="_")
colnames(calls_fp) = paste(colnames(calls_fp),"fp",sep="_")

calls_tp = calls_tp[order(sce$log10_total_counts,decreasing=FALSE),]
calls_tp = apply(calls_tp,2,cumsum)
df_tp    = data.frame(calls_tp)
df_tp$size = 1:nrow(df_tp)
df_tp = melt(df_tp, id.vars="size")
p_tp = ggplot(df_tp, aes(x=size,y=value,col=variable)) + geom_line() +
xlab("library size (rank)") + ylab("true positives (cummulative)") + xlim(0,5000) + ylim(0,90)


calls_fp = calls_fp[order(sce$log10_total_counts,decreasing=FALSE),]
calls_fp = apply(calls_fp,2,cumsum)
df_fp    = data.frame(calls_fp)
df_fp$size = 1:nrow(df_fp)
df_fp = melt(df_fp, id.vars="size")
p_fp = ggplot(df_fp, aes(x=size,y=value,col=variable)) + geom_line()+
xlab("library size (rank)") + ylab("false positives (cummulative)") + xlim(0,5000) + ylim(0,310)

pdf(file="./results/fig_dbl-scrublet-EDA.pdf",width=3*3*1.5,height=1*3*1.5)
plot_grid(p_tp,p_fp,nrow=1)
dev.off()
