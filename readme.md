# scds_manuscript

Instructions to create figures / tables from the [scds manuscript]().

### Step 1: Envrionment
We provide a docker container that contains all dependencies (python modules and R packages) necessary.

The following will set up the environment:

```{bash}
#- get docker image
$ docker pull kostkalab/scds:v2.0

#- run interactive session (bash), for instance:
$ docker run -it --name scds_tab_fig                                  \
                 --net=host                                           \
                 --env DISPLAY=$DISPLAY                               \
                 --volume $HOME/.Xauthority:/home/rstudio/.Xauthority \
                 kostkalab/scds:v2.0                                  \
                 /bin/bash
```

In the shell from above, get the code from github

```{bash}
$ cd
$ git clone https://github.com/kostkalab/scds_manuscript.git
$ cd scds_manuscript
```
### Step 2: Downloading and processing the data
The ```scds_manuscript/data``` folder contains directories for each of the four data
sets. Each directory contains two files, ```getRawData.sh``` and ```processRawData.R```.
```getRawData.sh``` downloads the data sets (or contains instructions on how to do so),
while ```processRawData.R``` can be run to process the input data and generate
```SingleCellExperiments``` that contain experimental doublet annotations. For each
data set, these annotations are stored in the ```colData``` slot.

For example, can  download the ```demuxlet``` data via

```{bash}
$ ./data/demu/getRawData.sh
```

The other data sets require human interaction (see the respective ```getRawData.sh``` files.

After downloading, the following generates the processed ```hgmm```
data:

```{bash}
$ [[ ! -d ./data/hgmm/proc ]] && mkdir ./data/hgmm/proc
$ R --vanilla < ./data/hgmm/processRawData.R
```
And analogous commands work for the other three data sets. Overall, this step generates four files containing processed data:

* ```./data/chcl/proc/sce_chcl.rds```
* ```./data/chpb/proc/sce_chpb.rds```
* ```./data/demu/proc/sce_demu.rds```
* ```./data/hgmm/proc/sce_hgmm.rds```

That all contain single cell experiments with experimental doublet annotations:

* ```colData(sce)$dbl_anno_lab``` Of type ```character``` for the experimental doublet annotation.
* ```colData(sce)$dbl_anno``` Of type ```logical``` for the experimental doublet annotation.

To process all data sets (after manual downloads where necessary) you can do:

```{bash}
for dat in {chcl,chpb,demu,hgmm}; do
  ./data/$dat/getRawData.sh
  [[ ! -d ./data/$dat/proc ]] && mkdir ./data/$dat/proc
  R --vanilla < ./data/$dat/processRawData.R
done
```


### Step 3: Computational doublet annotation

In this step we run ten doublet annotation methods on each of the four data sets.
For each method, we have a corresponding ```R``` script named ```run_[method].R``` in the ```./R``` subfolder.
So, ```./R/run_scrublet.R``` will be for running ```scrublet```, etc. There is also an ```apply_method.R``` script that takes two command line argumets - one providing the name of the method to use, the other the name of the data set to run the method on. So,

```{bash}
$ R --vanilla -f ./R/apply_method.R  --args --data=hgmm --method=scrublet
```

will run ```scrublet``` on the ```hgmm``` data, and

```{bash}
$ for method in {bcds,cxds,hybrid,scrublet,dblFinder,dblDetection,dblDecon,dblCells}; do
      for data in {chcl,chpb,demu,hgmm}; do
        R --vanilla -f ./R/apply_method.R  --args --data=$data --method=$method
    done
  done
```

will run four methods on all data sets (this will take a bit of time).

Results are written out in text files named according to ```./results/tmp/[dataset]_[method].txt```, which contain one line for each cell and three columns. The first
column contains the cell's barcode, the second the doublet score of the method, and the third column doublet calls (```TRUE```/```FALSE```, with ```TRUE``` indicating called doublets).
Higher scores indicate more confidence in annotating a doublet. For methods without calls (```dblCells```) or without scores (```dblDecon```) the respective entries are ```NA```.


Finally, the script ```./R/combine_results.R``` combines all doublet annotations for
each dataset and creates four ```SingleCellExperiments``` directly in the results directory:

```{bash}
$ R --vanilla < ./R/combine_results.R
```

* ```./results/sce_chcl.rds```
* ```./results/sce_chpb.rds```
* ```./results/sce_demu.rds```
* ```./results/sce_hgmm.rds```

Each data set now contains doublet annotations in the ```colData``` as columns with the following
names:

```{$method}_{score,call}```

where ```$method``` is one of ```{bcds, cxds, dblDecon, dblDetection, dblFinder, hybrid, scrublet, libsize, numfeat}```. The last two are baseline methods and have no associated calls. For example in ```R```:

```{R}
> sce.hgmm = readRDS("./results/sce_chpb.rds")
> sce.hgmm$hybrid_score  #- <<- hybrid scores
> sce.hgmm$scrublet_call #- <<- scrublet calls
```

This concludes doublet annotation.

### Step 4: Analyzing doublet annotation results
We analyze computational doublet annoations in four ```SingleCellExperiment``` objects, one for each data set (see Step 3 above).

#### Generate tables:

The ```mk_tabs.R``` script will generate latex tables in the ```results``` subdirectory.
They correspond to the tables in the manuscript as follows:

```{bash}
#- run R script to generate tables
$ R --vanilla < ./R/mk_tabs.R
$ cat ./results/tab_perf_avrg.tex
```

* ```tab_perf_avrg.tex``` Methods performance, averaged over data sets.
* ```tab_perf_cbnd.tex``` Methods performance for each data set in a sub-table.
* ```tab_dbldecon-comp.tex ``` Comparison of other methods to ```dblDecon``` for the ```ch_cell-lines``` data
* ```tab_hom-het.tex``` For the ```ch_cell-lines``` data, enrichment of hetertypic doublets in the true positive predictions for each method.
* ```tab_bcds-7.tex``` To compare heuristics for the number of training rounds in ```bcds```.
* ```tab_cxdsParams_ntop.tex``` To assess different values for the ```ntop``` parameter
* ```tab_cxdsParams_binThresh.tex``` To assess different values for the ```binThresh``` parameter.

#### Generate figures:

Figures are generated by the ```mk_figs.R``` script, analogously to the tables above:

```{bashr}
$ R --vanilla < ./R/mk_figs.R
```

 The figures generated are also written to the ```./results/``` directory. The following list of figures is produced:

* ```fig_cxds_[dataset].jpg```
* ```fig_perf_strat.pdf```
* ```fig_perf_strat_supp.pdf```
* ```fig_comp.pdf```
* ```fig_size-strat_[method].pdf```
* ```fig_tpfn_[method].jpg```
* ```fig_fdr-nfp.pdf```
* ```fig_fdr-fdr.pdf```

#### Assess running time:

To assess running time we limited resources with the following command:

```{bash}
$ cd /home/rstudio/scds_manuscript
$ taskset --cpu-list 2,3,19,20 R --vanilla < ./R/mk_running-times_python.R
$ cat ./results/tab_running-times.tex
```
#### Figure about robustness of performance assessment with respect to resampling of barcodes

Running the following creates the two resampling figures in the ```./results``` subdirectory.

```{bash}
$ nohup R --vanilla < ./R/mk_fig_resampling.R & #- takes a while
```
