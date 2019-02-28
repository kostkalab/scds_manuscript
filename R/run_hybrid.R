library(scds)
run_hybrid <- function(sce){
	sce = cxds_bcds_hybrid(sce)
	return(sce)
}


