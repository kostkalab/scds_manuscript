

#- fisher test for het enrichment in TPs:
getFET <- function(nme,sce){
#===========================

  hets = grepl("het",sce$type)
  pos  = sce$dbl_anno
  frac = mean(pos)
  tp   = pos & (binFu(colData(sce)[,nme],frac)==1)
  tab  = table(hets[pos],tp[pos])
  return(fisher.test(tab))
}

mk_hh_tab <- function(sce){
#==========================

    nms = c(  "cxds_score", "bcds_score", "hybrid_score",
              "scrublet_score", "dblFinder_score","dblDetection_score",
              "libsize_score", "numfeat_score", "dblCells_score")

    rlist = lapply(nms, getFET,sce)

    ors = unlist(lapply(rlist,function(x) x$estimate))
    pvs = unlist(lapply(rlist,function(x) x$p.value))
    tab = rbind(ors,pvs)
    colnames(tab) = nms
    tab = tab[,order(tab[2,],decreasing=FALSE)]
    rownames(tab)  = c("odds_ratio","p_value")
    tf  = tab
    tf[1,] = format(round(tab[1,],2),digits=2)
    tf[2,] = format(tab[2,],scientific=TRUE,digits=3)
    return(tf)
}

#- end
