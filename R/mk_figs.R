library(Matrix)
library(SingleCellExperiment)
library(scds)

#- READ IN DATA
#==============

sce.chcl = readRDS("./results/sce_chcl.rds")
sce.chpb = readRDS("./results/sce_chpb.rds")
sce.demu = readRDS("./results/sce_demu.rds")
sce.hgmm = readRDS("./results/sce_hgmm.rds")

#- FIGURE  (cxds gene pairs)
#============================
# fig_cxds

source("./R/fu_fig_cxds.R")

sce.chcl = addTSNE(sce.chcl)
sce.chpb = addTSNE(sce.chpb)
sce.demu = addTSNE(sce.demu)
sce.hgmm = addTSNE(sce.hgmm)

#- re-run cxds, because we lost the retRes info
#  by just saving the scores before

sce.chcl = cxds(sce.chcl, verb=TRUE, retRes=TRUE)
sce.chpb = cxds(sce.chpb, verb=TRUE, retRes=TRUE)
sce.demu = cxds(sce.demu, verb=TRUE, retRes=TRUE)
sce.hgmm = cxds(sce.hgmm, verb=TRUE, retRes=TRUE)

rowData(sce.demu)$cxds_symbol = tolower(rowData(sce.demu)$Symbol)
sce.demu = plotCXDSpairs(sce.demu,"./results/fig_cxds_demu.jpg")
rowData(sce.chpb)$cxds_symbol = tolower(rowData(sce.chpb)$gene)
sce.chpb = plotCXDSpairs(sce.chpb,"./results/fig_cxds_chpb.jpg")
rowData(sce.chcl)$cxds_symbol = tolower(rowData(sce.chcl)$gene)
sce.chcl = plotCXDSpairs(sce.chcl,"./results/fig_cxds_chcl.jpg")
rowData(sce.hgmm)$cxds_symbol = sub("hg19_","",tolower(rowData(sce.hgmm)$Symbol))
sce.hgmm = plotCXDSpairs(sce.hgmm,"./results/fig_cxds_hgmm.jpg")

#- FIGURE  (performance stratified by library size)
#===================================================
# fig_perf-strat

source("./R/fu_tab_performance.R")
source("./R/fu_fig_perf-strat.R")

pl.roc = list(  p1=stratPerf(sce.hgmm), p2=stratPerf(sce.demu),
                p3=stratPerf(sce.chpb), p4=stratPerf(sce.chcl))
pl.prc = list(  p1=stratPerf(sce.hgmm,AUROC=F), p2=stratPerf(sce.demu,AUROC=F),
                p3=stratPerf(sce.chpb,AUROC=F), p4=stratPerf(sce.chcl,AUROC=F))

filename="./results/fig_perf-strat.pdf"
pdf(file=filename,width=4*3*1.5,height=1*3)
plot_grid(plotlist=c(pl.roc[2],pl.prc[2],pl.roc[3],pl.prc[3]),ncol=4)
dev.off()

filename="./results/fig_perf-strat_supp.pdf"
pdf(file=filename,width=4*3*1.5,height=1*3)
plot_grid(plotlist=c(pl.roc[1],pl.prc[1],pl.roc[4],pl.prc[4]),ncol=4)
dev.off()

#- FIGURE  (upset plots)
#========================
# fig_comp

source("./R/fu_fig_comp.R") #- binFu already here.

pdf("/dev/null")
upSetFu(sce.hgmm)
grid.edit('arrange',name='arrange1')
p1 = grid.grab()
upSetFu(sce.demu)
grid.edit('arrange',name='arrange2')
p2 = grid.grab()
upSetFu(sce.chpb)
grid.edit('arrange',name='arrange3')
p3 = grid.grab()
upSetFu(sce.chcl)
grid.edit('arrange',name='arrange4')
p4 = grid.grab()
dev.off()
pdf(file="./results/fig_comp.pdf",width=2*15,height=5)
grid.arrange(p1,p2,p3,p4,nrow=1)
dev.off()

#- FIGURE  (density plots TP/FN/...)
#====================================
#- fig_tpfn

source("./R/fu_fig_tpfn.R")

makeTPFNplot(sce.chcl,"./results/fig_tpfn_chcl.jpg")
makeTPFNplot(sce.chpb,"./results/fig_tpfn_chpb.jpg")
makeTPFNplot(sce.demu,"./results/fig_tpfn_demu.jpg")
makeTPFNplot(sce.hgmm,"./results/fig_tpfn_hgmm.jpg")



#- FIGURE  library size stratified by TP/FP etc
#===============================================
#- fig_size-strat

source("./R/fu_fig_size-strat.R")

plotLibSzeStrat(sce.demu,filename="./results/fig_size-strat_demu.pdf")
plotLibSzeStrat(sce.chcl,filename="./results/fig_size-strat_chcl.pdf")
plotLibSzeStrat(sce.chpb,filename="./results/fig_size-strat_chpb.pdf")
plotLibSzeStrat(sce.hgmm,filename="./results/fig_size-strat_hgmm.pdf")

#- cummulative FPs / FDR by doublet score rank
#======================================

source("./R/fu_fig_FDR.R")

df = dfu(sce.chcl,"chcl")
df = rbind(df, dfu(sce.chpb,"chpb"))
df = rbind(df, dfu(sce.demu,"demu"))
df = rbind(df, dfu(sce.hgmm, "hgmm"))

theme_set(theme_gray())
p1 = ggplot(df,aes(x=rank,y=nfp,col=method)) + geom_line() + facet_wrap(~data, ncol=4,scales="free") +
xlab("doublet score rank") + ylab("# false positives")

pdf(file="./results/fig_fdr-nfp.pdf",width=4*3*1.5,height=1*3*1.5)
p1 %>% print
dev.off()

df = dfu(sce.chcl,"chcl",short=FALSE)
df = rbind(df, dfu(sce.chpb,"chpb",short=FALSE))
df = rbind(df, dfu(sce.demu,"demu",short=FALSE))
df = rbind(df, dfu(sce.hgmm, "hgmm",short=FALSE))

p2 = ggplot(df,aes(x=ftp,y=fdr,col=method)) + geom_line() + facet_wrap(~data, ncol=4) +
xlab("fraction of doublets recoverd") + ylab("false discorvery rate") +
xlim(0.1,1)

pdf(file="./results/fig_fdr-fdr.pdf",width=4*3*1.5,height=1*3*1.5)
p2 %>% print
dev.off()








#- end
