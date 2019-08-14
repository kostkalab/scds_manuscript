
#- we lmit resources either via docker or the following:
#- $taskset --cpu-list 2,3,7,9 /PATH/TO/R

library(SingleCellExperiment)
library(R.utils)
library(scds)
library(reticulate)
reticulate::use_python("/opt/user/rstudio/miniconda/bin/python", required = FALSE)
reticulate::use_condaenv("base",conda="/opt/user/rstudio/miniconda/bin/conda")
#reticulate::use_condaenv("R_reticulate_doublets")
sp = import("scipy.sparse")

source("./R/run_dblCells.R")
source("./R/run_dblFinder.R")

#================
#- PYTHON TIMINGS
#================
#
#- timing is done in python scripts, omitting data handling
#  here we just run them via system

PY_DD_SCRIPT = "./py/time_dblDetection.py" #- returns dd running time
PY_SC_SCRIPT = "./py/time_scrublet.py"     #- returns sc running time

time_py_meth <- function(cnts,nrep=3,timeout=180){
#=================================================

  #- write out data for python
  #---------------------------
  tfle = tempfile()
  message(tfle)
  sp$save_npz(tfle,cnts)
  npz_file = paste0(tfle,".npz")
  
  #- annotate doublets
  #-------------------
  cmd.dd = paste("python", PY_DD_SCRIPT, npz_file, "2>/dev/null | tail -n 1")
  cmd.sc = paste("python", PY_SC_SCRIPT, npz_file, "2>/dev/null | tail -n 1")

  out.sc = replicate(nrep,system(cmd.sc,intern=TRUE,timeout=timeout))
  out.sc = unlist(lapply(out.sc,function(x) if(length(x)==0) return(NA) else return(x)))
  sec.sc = as.integer(unlist(lapply(strsplit(out.sc,split=" \t "),"[[",2)))

  out.dd = replicate(nrep,system(cmd.dd,intern=TRUE,timeout=timeout))
  out.dd = unlist(lapply(out.dd,function(x) if(length(x)==0) return("foo \t NA") else return(x)))
  sec.dd = as.integer(unlist(lapply(strsplit(out.dd,split=" \t "),"[[",2)))
 
  #- clean up and return running time
  # (wall-clock, [sec, rounded to nearest integer], no data transfer (see the python scripts))
  file.remove(npz_file)
  res = cbind(sec.sc,sec.dd) ; colnames(res) = c("scrublet","dblDetection")
  
  return(res)
}

#===========
#- R TIMINGS
#===========

time_r_meth <- function(sce,nrep=3,timeout=180){
#===============================================
  
  #- NOTE: will not work without quoting expression/promise
  time_qexp = function(qexp){ 
                      replicate(nrep,withTimeout(round(system.time(eval(qexp),gcFirst=TRUE)["elapsed"]),
                                                 timeout = timeout, onTimeout = "silent" ))  
  }
  
  nafu          = function(lst) unlist(lapply(lst,function(x) if(is.null(x)) return(NA) else return(x)))
  sec.cxds      = nafu(time_qexp(quote(cxds(sce)         )))
  sec.bcds      = nafu(time_qexp(quote(bcds(sce)         )))
  sec.bcds7     = nafu(time_qexp(quote(bcds(sce,nmax=7)  )))
  sec.dblCells  = nafu(time_qexp(quote(run_dblCells(sce) ))) #somehow timeout does not work here properly...
  sec.dblFinder = nafu(time_qexp(quote(run_dblFinder(sce))))
    
  res = cbind(sec.cxds,sec.bcds,sec.bcds7,sec.dblCells,sec.dblFinder)
  colnames(res) = c("cxds","bcds","bcds7","dblCells","dblFinder")
  return(res)
}

#=====================
#- PUTTING IT TOGETHER
#=====================

mk_times <- function(sce,ncells,verb=TRUE){
#==========================================
  
  NREP   = 3
  TMEOUT = 60*30 #- half hour timeout
  ss     = sample(1:ncol(sce),ncells)
  tmp    = sce[,ss]
  tmp    = tmp[Matrix::rowSums(counts(tmp)>0)>0,]
  tmes_p = time_py_meth(counts(tmp),nrep=NREP,timeout=TMEOUT)
  tmes_r = time_r_meth(tmp,         nrep=NREP,timeout=TMEOUT)
  tmes   = apply(cbind(tmes_p,tmes_r),2,median)

  if(verb) message(ncells,":",tmes)
  #- wall-clock time, median over NREP replicates, TMEOUT sec max, rounded to nearest sec
  return(tmes) 
}

#===========
#- LOAD DATA
#===========

sce = readRDS("./data/demu/proc/sce_demu.rds")

#======================
#- ASSESS RUNNING TIMES
#======================

set.seed(124235)
sec.1k  = mk_times(sce, 1000)
sec.2k  = mk_times(sce, 2000)
sec.4k  = mk_times(sce, 4000)
sec.8k  = mk_times(sce, 8000)
sec.12k = mk_times(sce,12000)

#============
#- MAKE TABLE
#============

tab           = cbind(sec.1k,sec.2k,sec.4k,sec.8k,sec.12k)
colnames(tab) = c("1k","2k","4k","8k","12k")

tmp <- kable(tab[order(rowMeans(tab)),],"latex",booktabs=T) %>%
          kable_styling()
if(! dir.exists("./results")) dir.create("./results")
sink("./results/tab_running-times_python.tex")
print(tmp)
sink()

#- end
