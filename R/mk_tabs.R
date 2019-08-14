
library(Matrix)
library(PRROC)
library(pROC)
library(SingleCellExperiment)
library(knitr)
library(kableExtra)
library(scds)

#- READ IN DATA
#==============

sce.chcl = readRDS("./results/sce_chcl.rds")
sce.chpb = readRDS("./results/sce_chpb.rds")
sce.demu = readRDS("./results/sce_demu.rds")
sce.hgmm = readRDS("./results/sce_hgmm.rds")

#- PERFORMANCE METRICS ACROSS FOUR DATA SETS
#===========================================

#- functions to generate table
source("./R/fu_tab_performance.R")
#- files to write tables to, function for
pfix = "./results/tab_perf"
mk_tab <- function(sce, sfix){
  tab =  mk_perf_tab(sce)
  sink(paste(pfix,"_", sfix, ".tex",sep=""))
  print(kable(format(tab,digits=2),"latex"))
  sink()
  return(tab)
}
#- make tables
tab.chcl = mk_tab(sce.chcl,"chcl")
tab.chpb = mk_tab(sce.chpb,"chpb")
tab.demu = mk_tab(sce.demu,"demu")
tab.hgmm = mk_tab(sce.hgmm,"hgmm")
rm(mk_tab)

chcl_hl = paste("chcl: (", sum(sce.chcl$dbl_anno), " annotated doublets)", sep="")
chpb_hl = paste("chpb: (", sum(sce.chpb$dbl_anno), " annotated doublets)", sep="")
demu_hl = paste("demu: (", sum(sce.demu$dbl_anno), " annotated doublets)", sep="")
hgmm_hl = paste("hgmm: (", sum(sce.hgmm$dbl_anno), " annotated doublets)", sep="")

#- make one big table
sink("./results/tab_perf_cbnd.tex")

rbind(tab.chcl,tab.chpb,tab.demu,tab.hgmm)                %>%
round(2)                                                  %>%
kable(format="latex", booktabs=TRUE,linesep="",align="r") %>%
kable_styling()                                           %>%
group_rows(chcl_hl,1,9)                                   %>%
group_rows(chpb_hl,10,18)                                 %>%
group_rows(demu_hl,19,27)                                 %>%
group_rows(hgmm_hl,28,36)
sink()

#- average performance table:
rn  = rownames(tab.chcl)
tab = tab.chcl[rn,] + tab.chpb[rn,] + tab.demu[rn,] + tab.hgmm[rn,]
tab = round(tab/4,2)
tab = tab[order(tab[,1],tab[,8]),]
sink("./results/tab_perf_avrg.tex")
kable(round(tab,digits=2),format="latex", booktabs=TRUE, linesep = "")
sink()

#- COMPARISON WITH DOUBLET-DECON
#===============================

source("./R/fu_tab_dbldecon-comp.R")

dd.tab.chcl     = mk_dd_tab(sce.chcl)
dd.tab.chpb     = mk_dd_tab(sce.chpb)
dd.tab.demu     = mk_dd_tab(sce.demu)

tab = round(cbind(dd.tab.chcl,dd.tab.chpb,dd.tab.demu),2)
oo  = order(rowMeans(tab[,c(7,14,21)]))

sink("./results/tab_dbldecon-comp.tex")
tmp <- kable(tab[oo,c(1,5,6,7,8,12,13,14,15,19,20,21)],"latex",booktabs=T, linesep = "") %>%
         kable_styling() %>%
         add_header_above(c( " " = 1, "ch_cell-lines" = 4, "ch_pbmc" = 4, "demu"=4))
print(tmp)
sink()

#- HOM vs HET DOUBLETS IN SCE.CHCL
#=================================

#- better annotate hom vs het doublets:
#--------------------------------------
tmp  = as.character(sce.chcl$hto_classification)
tmp  = strsplit(tmp,split="_")
afu  = function(x) { if(length(x)<4)   return(x[1])
                     if(x[3] == x[1] ) return(paste(c("hom",sort(c(x[1],x[3]))),collapse="."))
                     return(paste(c("het",sort(c(x[1],x[3]))),collapse="."))}
mc   = unlist(lapply(tmp,afu))
mc   = tolower(mc)
sce.chcl$type = mc  # <<- new annotations for het vs. hom

source("./R/fu_tab_hom-het.R")
tab = mk_hh_tab(sce.chcl)
sink("./results/tab_hom-het.tex")
print(kable(tab,"latex", booktabs=TRUE))
sink()


#- COMPARISON cxds, bcds, bcds_7
#===============================

#- run bcds with nmax=7
sce.chcl.7 = bcds(sce.chcl, verb=TRUE, nmax=7)
sce.chpb.7 = bcds(sce.chpb, verb=TRUE, nmax=7)
sce.demu.7 = bcds(sce.demu, verb=TRUE, nmax=7)
sce.hgmm.7 = bcds(sce.hgmm, verb=TRUE, nmax=7)

ffu <- function(sce.7,sce){
  r1 = scorePred(sce$dbl_anno,sce$bcds_score)
  r2 = scorePred(sce$dbl_anno,sce$cxds_score)
  r3 = scorePred(sce.7$dbl_anno,sce.7$bcds_score)

  res = rbind(r1,r2,r3)
  rownames(res) = c("bcds","cxds","bcds_7")
  return(res)
}

t1  = ffu(sce.chcl.7,sce.chcl)
t2  = ffu(sce.chpb.7,sce.chpb)
t3  = ffu(sce.demu.7,sce.demu)
t4  = ffu(sce.hgmm.7,sce.hgmm)
tab = round(rbind(t1,t2,t3,t4),2)

sink("./results/tab_bcds-7.tex")
print(kable(rbind(tab)[,c(1:4,8)],"latex", booktabs=TRUE, linesep="") %>%
            kable_styling()                               %>%
            group_rows("ch_cell-lines",1,3)               %>%
            group_rows("ch_pbmc",4,6)                     %>%
            group_rows("demuxlet",7,9)                    %>%
            group_rows("hg-mm",10,12))
sink()

#- HERE we count NA's per method
#===============================

countNA <- function(sce){

  res = c(  sum(is.na(sce$cxds_score)),
            sum(is.na(sce$bcds_score)),
            sum(is.na(sce$hybrid_score)),
            sum(is.na(sce$scrublet_score)),
            sum(is.na(sce$dblFinder_score)),
            sum(is.na(sce$dblCells_score)),
            sum(is.na(sce$dblDetection_score)),
            sum(is.na(sce$libsize_score)),
            sum(is.na(sce$numfeat_score)))
}

rbind(  countNA(sce.chcl),
        countNA(sce.chpb),
        countNA(sce.demu),
        countNA(sce.hgmm))


#===========================================
#- Do cxds parameters influence performance?
#===========================================

source("./wrk_tab_cxds-params.R")






#- end
