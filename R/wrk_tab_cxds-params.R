#- START TBR
library(reshape2)
library(dplyr)
library(kableExtra)

sce.chcl = readRDS("./results/sce_chcl.rds")
sce.chpb = readRDS("./results/sce_chpb.rds")
sce.demu = readRDS("./results/sce_demu.rds")
sce.hgmm = readRDS("./results/sce_hgmm.rds")

#- need performance measures
source("./R/fu_tab_performance.R")

set.seed(23124)

#- data sets and parameter ranges
dats           <- c("chcl", "chpb", "demu", "hgmm")
binThresh_vals <- seq(0, 5, by = 1)
ntop_vals      <- seq(200, 1000, by = 100)

RES.NTOP      = NULL
RES.BINTHRESH = NULL

for(dat in dats){

    #- read in data set
    data_file     = paste("./data/",dat,"/proc/sce_",dat,".rds",sep="")
    sce           = readRDS(data_file)

    #- predict with thresholds and evaluate (ntop=500)
    afu.bt <- function(binThresh,sce) {
        sce = cxds(sce,verb=TRUE,binThresh=binThresh)
        res = scorePred(resp=sce$dbl_anno,pred=sce$cxds_score)
        return(res)
    }
    res.binThresh           = sapply(binThresh_vals,afu.bt,sce)
    colnames(res.binThresh) = binThresh_vals
    res.binThresh           = melt(res.binThresh)
    colnames(res.binThresh) = c("metric","binThresh","value")
    res.binThresh$data      = dat
    res.binThresh$ntop      = 500

    #- predict with different ntops and evaluate (binThresh=0)
    afu.nt <- function(ntop,sce) {
        sce = cxds(sce,verb=TRUE,ntop=ntop)
        res = scorePred(resp=sce$dbl_anno,pred=sce$cxds_score)
        return(res)
    }
    res.ntop           = sapply(ntop_vals,afu.nt,sce)
    colnames(res.ntop) = ntop_vals
    res.ntop           = melt(res.ntop)
    colnames(res.ntop) = c("metric","ntop","value")
    res.ntop$data      = dat
    res.ntop$binThresh = 500

    #- keep track across data sets
    RES.NTOP      = rbind(RES.NTOP, res.ntop)
    RES.BINTHRESH = rbind(RES.BINTHRESH, res.binThresh)
}

#- organize data
tab.nt      = data.frame(round(t(acast(RES.NTOP,metric~data+ntop)),3))
tab.nt$ntop = unlist(lapply(strsplit(rownames(tab.nt),split="_"),"[[",2))
tab.nt$data = unlist(lapply(strsplit(rownames(tab.nt),split="_"),"[[",1))
tab.bt           = data.frame(round(t(acast(RES.BINTHRESH,metric~data+binThresh)),3))
tab.bt$binThresh = unlist(lapply(strsplit(rownames(tab.bt),split="_"),"[[",2))
tab.bt$data      = unlist(lapply(strsplit(rownames(tab.bt),split="_"),"[[",1))

tab.nt.s <- tab.nt %>% group_by(data)                    %>% select(data,ntop,AUC,AUPRC,pAUC950,pAUC990,rec05,prec05) %>%
arrange(data,desc(AUC),desc(AUPRC))                      %>%
kable(format="latex",booktabs=TRUE)                      %>%
kable_styling()                                          %>%
group_rows("chcl",1,9)                                   %>%
group_rows("chpb",10,18)                                 %>%
group_rows("demu",19,27)                                 %>%
group_rows("hgmm",28,36)

tab.bt.s <- tab.bt %>% group_by(data)                         %>% select(data,binThresh,AUC,AUPRC,pAUC950,pAUC990,rec05,prec05) %>%
arrange(data, desc(AUC),desc(AUPRC))                          %>%
kable(format="latex",booktabs=TRUE)                           %>%
kable_styling()                                               %>%
group_rows("chcl",1,6)                                        %>%
group_rows("chpb",7,12)                                       %>%
group_rows("demu",13,18)                                      %>%
group_rows("hgmm",19,24)

sink("./results/tab_cxdsParams_ntop.tex")
print(tab.nt.s)
sink()

sink("./results/tab_cxdsParams_binThresh.tex")
print(tab.bt.s)
sink()
