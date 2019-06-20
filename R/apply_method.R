library(SingleCellExperiment)

#-- parse arguments
ca  = commandArgs(trailingOnly = TRUE)
dat = gsub("--data=",'',ca[grepl("--data",ca)])
met = gsub("--method=",'',ca[grepl("--method",ca)])

#- files and functions
method_script = paste("./R/run_",met,".R",sep="")
method_name   = paste("run_",met,sep="")
data_file     = paste("./data/",dat,"/proc/sce_",dat,".rds",sep="")
out_file      = paste("./results/tmp/",dat,"_",met,".txt",sep="")

#- don't overwrite
if(file.exists(out_file)) {
	cat("\n Output file already exists.\n exiting.\n")
	q(save="no")
}

#- run doublet annotation
source(method_script)
sce = readRDS(data_file)
sce = do.call(method_name,list(sce=sce))

#- save annotation scores and calls
if(!dir.exists("./results")) dir.create("./results")
if(!dir.exists("./results/tmp")) dir.create("./results/tmp")
write.table(cbind( colnames(sce),
                   colData(sce)[,paste(met,"score",sep="_")],
                   colData(sce)[,paste(met,"call",sep="_")]),
	          file=out_file,quote=FALSE,sep="\t",col.names=FALSE,row.names=FALSE)
