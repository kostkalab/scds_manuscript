
#- we lmit resources either via docker or the following:
#- $taskset --cpu-list 2,3 /PATH/TO/R

library(microbenchmark)
library(S4Vectors)
library(R.utils) #- evaluate with timeout
library(reshape2)
library(kableExtra)

source("./R/run_cxds.R")
source("./R/run_bcds.R")
source("./R/run_scrublet.R")
source("./R/run_dblDetection.R")
source("./R/run_dblFinder.R")
source("./R/run_dblCells.R")

set.seed(124235)

#- FIXME: scrublet reports faster times - wrapper issue?

sce = readRDS("./data/demu/proc/sce_demu.rds")

timeRun <- function(fu,sce,DF,DFind,name,...){
#=============================================

  #- run 3 times and take middle value; abort when running time > 5min
  #- still need tryCatch because withTimeout has error with reticulate
  tmp = tryCatch( withTimeout( { microbenchmark(a=fu(sce,...),times=3L,unit="s") },
                       timeout = 3*10*60, onTimeout = "silent" ),
                  error = function(e) return(NULL) )
  DF[DFind,1] = name
  if(is.null(tmp)){
    DF[DFind,2] = NA
  } else {
    DF[DFind,2] = median(tmp$time)/1000/1000/1000
  }
  DF[DFind,3] = ncol(sce)
  return(DF)
}


makeDF <- function(sce){
#=======================

  DF = DataFrame(method=character(),time=numeric(),ncells=integer())
  DF = timeRun(run_cxds,sce,DF,1,"cxds")
  DF = timeRun(bcds,sce,DF,2,"bcds_7",nmax=7)
  DF = timeRun(bcds,sce,DF,3,"bcds",nmax="tune")
  DF = timeRun(run_scrublet,sce,DF,4,"scrublet")
  DF = timeRun(run_dblDetection,sce,DF,5,"dblDetection")
  DF = timeRun(run_dblFinder,sce,DF,6,"dblFinder")
  DF = timeRun(run_dblCells,sce,DF,7,"dblCells")

  return(DF)
}

ss  = sample(1:ncol(sce),1000)
sces = sce[,ss]
DF.1k = makeDF(sces)

ss  = sample(1:ncol(sce),2000)
sces = sce[,ss]
DF.2k = makeDF(sces)

ss  = sample(1:ncol(sce),4000)
sces = sce[,ss]
DF.4k = makeDF(sces)

ss  = sample(1:ncol(sce),8000)
sces = sce[,ss]
DF.8k = makeDF(sces)

ss  = sample(1:ncol(sce),12000)
sces = sce[,ss]
DF.12k = makeDF(sces)

df  = as.data.frame(rbind(DF.1k,DF.2k,DF.4k, DF.8k,DF.12k))
tab = round(acast(df,method ~ ncells,value.var="time"),1)

tmp <- kable(tab[order(rowMeans(tab)),],"latex",booktabs=T) %>%
          kable_styling()
if(! dir.exists("./results")) dir.create("./results")
sink("./results/tab_running-times.tex")
print(tmp)
sink()

#- end
