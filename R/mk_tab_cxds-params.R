#- START TBR
source("./R/wrk_tab_cxds-params.R")

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
