
library(UpSetR)
library(grid)
library(gridExtra)

source("./R/fu_fig_cxds.R") #- binfu

upSetFu <- function(sce){
#========================

dat               = data.frame(call = sce$dbl_anno)
dat$annotation    = as.integer(dat$call)
frc               = mean(dat$call)
dat$cxds          = binFu(sce$cxds_score,frc)
dat$scrublet      = binFu(sce$scrublet_score,frc)
dat$bcds          = binFu(sce$bcds_score,frc)
dat$dblFinder     = binFu(sce$dblFinder_score,frc)
dat$dblDetection  = binFu(sce$dblDetection_score,frc)
dat$hybrid        = binFu(sce$hybrid_score,frc)
dat$libsize       = binFu(sce$libsize_score,frc)
dat$features      = binFu(sce$numfeat_score,frc)
dat$dblCells      = binFu(sce$dblCells_score,frc)

hasAtt <- function(row,att) newData <- (row[att] >0)

p = upset(dat,nsets=10,nintersects=20,show.numbers=FALSE,
  main.bar.color=rgb(0,0,0,1/3), matrix.color=rgb(0,0,0,1), set_size.show=FALSE,
  order.by=c("freq"), text.scale=1.4,mb.ratio=c(0.6,0.4),matrix.dot.alpha=0,
 queries = list(list(query = hasAtt, params = list("call"), color=rgb(0,0,0,1), active=TRUE))
)

#- remove bars manually

return(p)
}
