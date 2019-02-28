
library(SingleCellExperiment)
library(Matrix)

dats = c("chcl","chpb","demu","hgmm")
mets = c("bcds","cxds","dblCells","dblDecon","dblDetection","dblFinder","hybrid","scrublet")

for(da in dats){
	cat("\n Combining data set: ", da, "\n")
	sce_file = paste("./data/",da,"/proc/sce_",da,".rds",sep="")
	out_file = paste("./results/sce_",da,".rds",sep="")
	#- don't overwrite
	if(file.exists(out_file)){
		cat("\n Output file already exists.\n exiting.\n")
		# cat(outfile)
		cat("\n")
        	q(save="no")
	}
	sce      = readRDS(sce_file)
	for(me in mets){
		cat("    Combining method: ", me, "\n")
		anno_file = paste("./results/tmp/",da,"_",me,".txt",sep="")
		if(! file.exists(anno_file)){
			cat("    --> Missing annotation file! Skipping. \n")
			next
		}
		anno = read.table(anno_file,stringsAsFactors=FALSE)
		scrs = anno[,2] ; names(scrs) = anno[,1]
		if(! all(colnames(sce) %in% names(scrs)) ){
			cat("    --> Missing barcodes! Skipping. \n")
      next
  	}
    if(me != "dblDecon"){
		      cn = paste(me,"_","score",sep="")
		      colData(sce)[,cn] = scrs[colnames(sce)]
    } else {
      cn = paste(me,"_","call",sep="")
      colData(sce)[,cn] = scrs[colnames(sce)]
    }
	}

  #- add baseline methods
  sce$libsize_score = Matrix::colSums(counts(sce))
  sce$numfeat_score = Matrix::colSums(counts(sce)>0)
	saveRDS(sce,out_file)
}
