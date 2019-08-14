

library(SingleCellExperiment)
library(knitr)
library(kableExtra)
library(dplyr)
source("./R/fu_tab_calls.R")

sce.chcl = readRDS("./results/sce_chcl.rds")
sce.chpb = readRDS("./results/sce_chpb.rds")
sce.demu = readRDS("./results/sce_demu.rds")
sce.hgmm = readRDS("./results/sce_hgmm.rds")

tab.hgmm = mk_call_tab(sce.hgmm); tab.hgmm$data = "hgmm"
tab.demu = mk_call_tab(sce.demu); tab.demu$data = "demu"
tab.chpb = mk_call_tab(sce.chpb); tab.chpb$data = "chpb"
tab.chcl = mk_call_tab(sce.chcl); tab.chcl$data = "chcl"

tab.hgmm <- tab.hgmm %>% filter(method != "dblDecon")

rbind(tab.chcl,tab.chpb,tab.demu,tab.hgmm)        %>%
  filter(cut_type %in% c("built_in","balanced"))                    %>%
  mutate(f1=round(2*precision*recall/(precision+recall),2))         %>%
  group_by(data)                                                    %>%
  arrange(data,desc(f1),desc(balAcc),desc(precision),desc(recall), desc((num_dbl)))  %>%
  ungroup                                                           %>%
  select(num_dbl,method,f1,balAcc,precision,recall,FN)              %>%
  kable(format="latex",booktabs=TRUE)                               %>%
  kable_styling()                                                   %>%
  group_rows("chcl (889 doublets)",1,7)                             %>%
  group_rows("chpb (2,598 doublets)",8,14)                          %>%
  group_rows("demu (1,565 doublets)",15,21)                         %>%
  group_rows("hgmm (741 doublets)",22,27)


#==========================
#- false positive data sets
#==========================

#- NOTE: data not included in docker image.
#

sce.micm = readRDS("./sce_micm.rds")
sce.bmpr = readRDS("./sce_bmpr.rds")

runMeth <- function(method_name, sce){
  method_script = paste("./R/run_",method_name,".R",sep="")
  function_name = paste("run_",method_name,sep="")
  source(method_script)
  sce <- do.call(function_name,list(sce=sce))
  return(sce)
}

#- annotate
#==========

sce.bmpr <- runMeth("bcds",sce.bmpr)
sce.micm <- runMeth("bcds",sce.micm)

sce.bmpr <- runMeth("cxds",sce.bmpr)
sce.micm <- runMeth("cxds",sce.micm)

sce.bmpr <- runMeth("hybrid",sce.bmpr)
sce.micm <- runMeth("hybrid",sce.micm)

sce.bmpr <- runMeth("dblDetection",sce.bmpr)
sce.micm <- runMeth("dblDetection",sce.micm)

sce.bmpr <- runMeth("dblFinder",sce.bmpr)
sce.micm <- runMeth("dblFinder",sce.micm)

sce.bmpr <- runMeth("scrublet",sce.bmpr)
sce.micm <- runMeth("scrublet",sce.micm)

#- make table
#============
sce.bmpr$dbl_anno = FALSE
sce.micm$dbl_anno = FALSE
sce.bmpr$dblDecon_call = NA
sce.micm$dblDecon_call = NA

tab.bmpr = mk_call_tab(sce.bmpr); tab.bmpr$data = "bmpr"
tab.micm = mk_call_tab(sce.micm); tab.micm$data = "micm"

tab.bmpr <- tab.bmpr %>% filter(method != "dblDecon")
tab.micm <- tab.micm %>% filter(method != "dblDecon")

rbind(tab.bmpr,tab.micm)        %>%
  arrange(data,desc(num_dbl))   %>%
  ungroup                                                           %>%
  select(num_dbl,method)              %>%
  kable(format="latex",booktabs=TRUE)                               %>%
  kable_styling()                                                   %>%
  group_rows("bmpr (0 doublets)",1,6)                             %>%
  group_rows("micm (0 doublets)",7,12)





#- end
