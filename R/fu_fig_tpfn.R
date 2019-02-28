
makeTPFNplot <- function(sce, filename){
#======================================


    dat              = reducedDim(sce,"tsne")
    colnames(dat)    = c("x","y")
    dat              = as.data.frame(dat)
    dat$dblDetection = rank(sce$dblDetection_score)
    dat$dblFinder    = rank(sce$dblFinder_score)
    dat$hybrid       = rank(sce$hybrid_score)
    dat$bcds         = rank(sce$bcds_score)
    dat$cxds         = rank(sce$cxds_score)
    dat$scrublet     = rank(sce$scrublet_score)
    dat$dblCells     = rank(sce$dblCells_score)
    dat$libsize      = rank(sce$libsize_score)
    dat$features     = rank(sce$numfeat_score)
    dat$annotation   = as.integer(sce$dbl_anno)


      makeCol <- function(cn,dat,leOnly=FALSE){
      #------------------------------------
      resp =  as.integer(dat$annotation==1)
      pred =  binFu(dat[,cn],mean(resp))

          dat$TP = resp*pred
          dat$FP = (1-resp)*pred
          dat$FN = resp*(1-pred)

          p.sc = plotPoints(cn,dat,au=.5)
          p.sc = p.sc +  guides(alpha = guide_legend(override.aes = list(size=10),title="Rank")) +
                 theme(legend.position="right",legend.key.size=unit(2.5,"line"),
                 legend.text=element_text(size=35),legend.title=element_text(size=35))
          le.sc = get_legend(p.sc)
          p.sc = p.sc + theme(legend.position="none")
          x.range = ggplot_build(p.sc)$layout$panel_params[[1]]$x.range
          y.range = ggplot_build(p.sc)$layout$panel_params[[1]]$y.range

          p.tp = ggplot(dat[dat$TP==1,], aes(x=x,y=y)) +
                 stat_density_2d(geom = "tile", aes(fill = stat(density)) ,
                 contour = FALSE) + scale_fill_gradient(low= rgb(1,1,1),
                 high = "olivedrab") +
                 geom_point(aes(x=x,y=y),size=.1,color="black")+
                 theme(line=element_blank(), panel.background = element_blank(), axis.text.x=element_blank(),
                 axis.text.y=element_blank(),      axis.title.x=element_blank(),axis.title.y=element_blank()) +
                 guides(fill = guide_colourbar(label = FALSE, ticks=FALSE,)) +
                 theme(legend.position="right",legend.key.size=unit(2.5,"line"),
                 legend.text=element_text(size=35),legend.title=element_text(size=35))
          le.tp = get_legend(p.tp)
          p.tp = p.tp+theme(legend.position="none") + xlim(x.range[1],x.range[2]) +
                 ylim(y.range[1],y.range[2])


          p.fp = ggplot(dat[dat$FP==1,], aes(x=x,y=y)) +
                stat_density_2d(geom = "tile", aes(fill = stat(density)) ,
                contour = FALSE) + scale_fill_gradient(low= rgb(1,1,1),
                high = "firebrick") +
                geom_point(aes(x=x,y=y),size=.1,color="black")+
                theme(line=element_blank(), panel.background = element_blank(), axis.text.x=element_blank(),
                axis.text.y=element_blank(),      axis.title.x=element_blank(),axis.title.y=element_blank()) +
                guides(fill = guide_colourbar(label = FALSE, ticks=FALSE,)) +
                theme(legend.position="right",legend.key.size=unit(2.5,"line"),
                legend.text=element_text(size=35),legend.title=element_text(size=35))
          le.fp = get_legend(p.fp)
          p.fp= p.fp+theme(legend.position="none") + xlim(x.range[1],x.range[2]) +
                 ylim(y.range[1],y.range[2])

          p.fn = ggplot(dat[dat$FN==1,], aes(x=x,y=y)) +
               stat_density_2d(geom = "tile", aes(fill = stat(density)) ,
               contour = FALSE) + scale_fill_gradient(low= rgb(1,1,1),
               high = "steelblue") +
               geom_point(aes(x=x,y=y),size=.1,color="black")+
               theme(line=element_blank(), panel.background = element_blank(), axis.text.x=element_blank(),
               axis.text.y=element_blank(),      axis.title.x=element_blank(),axis.title.y=element_blank()) +
               guides(fill = guide_colourbar(label = FALSE, ticks=FALSE,))+
               theme(legend.position="right",legend.key.size=unit(2.5,"line"),
               legend.text=element_text(size=35),legend.title=element_text(size=35))
          le.fn = get_legend(p.fn)
          p.fn = p.fn+theme(legend.position="none") + xlim(x.range[1],x.range[2]) +
                 ylim(y.range[1],y.range[2])

         if(leOnly) return(list(le.sc,le.tp,le.fp,le.fn))

         return(list(p.sc,p.tp,p.fp,p.fn))
     }
     pdf("/dev/null")
     plst =  lapply(colnames(dat)[3:11],makeCol,dat)
     dev.off()
     jpeg(file=filename,width=4*300,height=10*300,pointsize=60)
     print(plot_grid(plotlist=c(plst[[1]],plst[[2]],plst[[3]],plst[[4]],plst[[5]],
       plst[[6]], plst[[7]],plst[[8]],plst[[9]], makeCol(colnames(dat)[9],dat,leOnly=T)),ncol=4))
     dev.off()
 }
