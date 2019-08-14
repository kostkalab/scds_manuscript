library(scds)
run_bcds <- function(sce){
	sce = bcds(sce,verb=TRUE,estNdbl=TRUE)
	return(sce)
}
