library(Matrix)
library(PRROC)
library(pROC)
library(SingleCellExperiment)
library(knitr)
library(kableExtra)
library(scds)
library(pamr)

#- READ IN DATA
#==============

sce.chcl = readRDS("./results/sce_chcl.rds")
sce.chpb = readRDS("./results/sce_chpb.rds")
sce.demu = readRDS("./results/sce_demu.rds")
sce.hgmm = readRDS("./results/sce_hgmm.rds")



#- Can we learn a combination of dbl scores?
#===========================================
library(pamr)
source("./R/fu_tab_combination-learning.R")

dat.hgmm = getDAT(sce.hgmm)
dat.demu = getDAT(sce.demu)
dat.chcl = getDAT(sce.chcl)
dat.chpb = getDAT(sce.chpb)

mod.chcl = getModel(dat.chcl$X,dat.chcl$y)
mod.chpb = getModel(dat.chpb$X,dat.chpb$y)
mod.demu = getModel(dat.demu$X,dat.demu$y)
mod.hgmm = getModel(dat.hgmm$X,dat.hgmm$y)

chcl.chcl = assessModelCV(dat.chcl$X,dat.chcl$y)
chpb.chpb = assessModelCV(dat.chpb$X,dat.chpb$y)
demu.demu = assessModelCV(dat.demu$X,dat.demu$y)
hgmm.hgmm = assessModelCV(dat.hgmm$X,dat.hgmm$y)

mods = list(mod.chcl,mod.chpb,mod.demu,mod.hgmm)
dats = list(dat.chcl,dat.chpb,dat.demu,dat.hgmm)

bst_single = c( assessBestFeature(dat.chcl),
                assessBestFeature(dat.chpb),
                assessBestFeature(dat.demu),
                assessBestFeature(dat.hgmm))


res           = apply(expand.grid(1:4,1:4),1,function(x) assessModel(mods[[x[1]]],dats[[x[2]]]$X,dats[[x[2]]]$y))
res           = matrix(res,ncol=4,byrow=TRUE)
res           = cbind(res,bst_single)
colnames(res) = paste("mod",c("chcl","chpb","demu","hgmm","bst_single"),sep="_")
rownames(res) = paste("dat",c("chcl","chpb","demu","hgmm"),sep="_")
diag(res)     = c(chcl.chcl,chpb.chpb,demu.demu,hgmm.hgmm)
kable(round(res,2),format="pandoc")
