
library(scran)
library(rsvd)
library(Rtsne)
library(cowplot)
library(ggplot2)


addTSNE <- function(sce){
#========================
  set.seed(32141)
  lc  = t(log1p(counts(sce)))
  lc  = lc/Matrix::rowMeans(lc)
  vr  = apply(lc,2,var)
  lc  = lc[,order(vr,decreasing=TRUE)[1:500]]
  pc  = rpca(lc,k=10)
  rdr = Rtsne(pc$x,verb=T)
  reducedDim(sce,"tsne") = rdr$Y

  return(sce)
}


binFu <- function(x,frc){
#========================

  numpos = round(length(x)*frc)
  xin    = order(x,decreasing=TRUE)[1:numpos]
  res    = as.integer(1:length(x) %in% xin)
  return(res)
}


plotCoex <- function(g1,g2,dat){
#===============================

    #- coexpression (average)
    dat$coex = (dat[,g1]+dat[,g2])/2
    dat$coex[dat[,g1]==0] = 0
    dat$coex[dat[,g2]==0] = 0
    p = ggplot(data=dat,aes_string(x="x",y="y",alpha="coex"))+ scale_alpha_continuous(range=c(0,1)) + geom_point(size=.4) +
    theme(line=element_blank(), panel.background = element_blank(), axis.text.x=element_blank(),
    axis.text.y=element_blank(),      axis.title.x=element_blank(),axis.title.y=element_blank()) +
    theme(legend.position="none")+ggtitle(paste("Coexpression")) +
    theme(plot.title = element_text(size=28,color=rgb(0,0,0,.7)))
    return(p)
}

plotPoints <- function(vname,dat,col=rgb(0,0,0,1),al=0.0,au=1,title=""){
#=======================================================================

    p = ggplot(data=dat,aes_string(x="x",y="y",alpha=vname))+ scale_alpha_continuous(range=c(al,au)) +
    geom_point(size=.5,color=col) +
    theme(line=element_blank(), panel.background = element_blank(), axis.text.x=element_blank(),
    axis.text.y=element_blank(),      axis.title.x=element_blank(),axis.title.y=element_blank()) +
    theme(legend.position="none")+ggtitle(title)+
    theme(plot.title = element_text(size=28,color=rgb(0,0,0,.7)))
    return(p)
}


plotCXDSpairs <- function(sce,filename){
#=======================================

      if(is.null(metadata(sce)$cxds_topPairs)) stop("Error: No gene pairs found\n")

      #- might need to run tsne
      if(is.null(reducedDim(sce,"tsne"))) stop("Error: tSNE projection found\n")

      #- Generate data frame used for plotting
      #---------------------------------------

      #- gene names:
      if(is.null(rowData(sce)$cxds_symbol)){
        if(is.null(rownames(sce))) {
          stop("Error: sce needs rownames. Please fix\n")
        } else {
          rowData(sce)$cxds_symbol = rownames(sce)
        }
      }
      rs        = rowData(sce)$cxds_symbol
      hb        = rowData(sce)$cxds_hvg_bool
      ho        = rowData(sce)$cxds_hvg_ordr[hb]
      hgs       = rs[ho]
      gP        = apply( metadata(sce)$cxds_topPairs, 1, function(x) hgs[x])
      gs        = as.vector(gP)[1:50]
      dat       = data.frame(reducedDim(sce,"tsne"))
      names(ho) = hgs
      dat       = cbind(dat,as.matrix(t(log1p(counts(sce)[ho[unique(gs)],]))))
      gs        = tolower(gsub("-",".",gs))
      gs        = tolower(gsub("\\|",".",gs))

      colnames(dat)    = c("x","y",unique(gs))
      dat$doublets     = sce$dbl_anno
      dat$predictions  = binFu(sce$cxds_score, mean(sce$dbl_anno) )
      dat$Doublets     = as.numeric(dat$doublets)
      dat$Predictions  = as.numeric(dat$predictions)
      dat$Cells        = 1

      #make list of plots
      makeRow = function(i){
      list(plotPoints(gs[i],dat,title=gs[i]),plotPoints(gs[i+1],dat,title=gs[i+1]),plotCoex(gs[i],gs[i+1],dat))
      }
      plst = list( plotPoints("Cells",      dat,col=rgb(.2,.2,.7),au=.3,title="Cells"),
                   plotPoints("Doublets",   dat,col=rgb(.7,.2,.2),au=.3,title="Doublets"),
                   plotPoints("Predictions",dat,col=rgb(.7,.2,.2),au=.3,title="Predictions"))
      plst = c(plst, makeRow(1),makeRow(3),makeRow(5),makeRow(7),makeRow(9))

      jpeg(file=filename,width=3*300,height=6*300,pointsize=60)
      print(plot_grid(plotlist=plst,ncol=3))
      dev.off()
      return(sce)
  }
