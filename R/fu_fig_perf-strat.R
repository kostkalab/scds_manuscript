library(gtools)
library(stringr)
library(ggplot2)
library(reshape2)


stratPerf <- function(sce, leg=TRUE, AUROC=TRUE){
#================================================

  bn = quantcut( sce$libsize_score, 10)
  if(AUROC){
    pind = 1
  } else {
    pind = 8
  }

  ff  = function(sce) mk_perf_tab(sce,ord=FALSE)[,pind]
  acs = tapply(1:length(bn),bn,function(x) ff(sce[,x]))
  mat = matrix(unlist(acs),nrow=9)

  colnames(mat) = str_pad(1:10,2,pad="0")
  #colnames(mat) = paste("quantile_",1:10,sep="")
  rownames(mat) = c("cxds","bcds","hybrid","libsize","features",
                    "scrublet","dblDetection","dblFinder","dblCells")

  #- plotting
  if(AUROC){
    dfm = melt(mat,measured.vars=1:10,value.name="AUROC")
    vs  = "AUROC"
  } else {
    dfm = melt(mat,measured.vars=1:10,value.name="AUPRC")
    vs  = "AUPRC"
  }
  colnames(dfm)[1:2] = c("method","quantile")
  dfm$quantile = str_pad(dfm$quantile,2,pad="0")
  dfm$method = factor(dfm$method, levels = levels(dfm$method)[order(rowMeans(mat))])
  p = ggplot(dfm, aes(quantile, method)) + geom_tile(aes_string(fill = vs),
       colour = "white") + scale_fill_gradient(low = "white", high = "steelblue")+# limits=c(0.4,1))+
       theme_bw() + labs(x = "",y = "")  +
       theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=12), axis.text.y = element_text(size=12)) +
       xlab("") + ylab("")
  if(leg==FALSE) p = p+theme(legend.position="none")
  return(p)
}
