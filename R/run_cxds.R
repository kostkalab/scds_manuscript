library(scds)
run_cxds <- function(sce){
	sce = cxds(sce,verb=TRUE,retRes=TRUE)
	return(sce)
}
