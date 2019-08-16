# Readme for misc/ folder

The folder "negCtrl_data/" contains code to retrieve 2 "negative" control datasets,
where doublets have been visually inspected and excluded. These are needed for some additional analyses provided in file: R/calls.R.

* ```micm```:  This dataset contains cells from mouse pre-implantation inner cell mass (ICM), epiblast and extra-embryonic endoderm across various stages of embryonic development (https://doi.org/10.1016/j.celrep.2017.07.009), where putative doublets were identified via visual inspection during cell picking. Data files for gene expression counts were downloaded from GEO \href{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE100597}{(GSE100597)}. Excluding putative doublets yielded a set of 24,483 genes and 588 singlet cells.


* ```bmpr```: This dataset contains hematopoietic bone marrow progenitor cells from mouse where doublets have been identified and excluded via cell capture imaging (https://doi.org/10.1038/nature19348). It contains 382 cells, consisting of CMP, GMP, LK CD34 and LSK cells. This data is also mentioned in the doubletDecon preprint (https://doi.org/10.1101/364810). To retrieve raw gene expression counts, we downloaded fastq files from GEO \href{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70245}{(GSE70245)}. Reads were aligned to the mouse GRCm38.p6 genome using HISAT2 (v2.1.0) and gene expression counts retrieved using the ```featureCounts``` function from the R package ```Rsubread``` (v1.32.4) with the Gencode vM21 gtf file. This yielded data comprising expression counts for 26,829 genes across 382 singlet cells.
