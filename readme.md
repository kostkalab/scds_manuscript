# scds_manuscript

Instructions to create figures / tables from the [scds manuscript]().

### Step 1: Envrionment
We provide a docker container that contains all dependencies (python modules and R packages) necessary.

The following will set up the environment:

```{bash}
#- get docker image
$ docker pull kostkalab/scds:v1.0

#- run interactive session (bash), for instance:
$ docker run -it --name scds_tab_fig                                  \
                 --net=host                                           \
                 --env DISPLAY=$DISPLAY                               \
                 --volume $HOME/.Xauthority:/home/rstudio/.Xauthority \
                 kostkalab/scds:v1.0                                  \
                 /bin/bash
```

In the shell from above, get the code from github

```{bash}
$ cd
$ git clone https://dkostka@bitbucket.org/dkostka/scds_manuscript.git
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
$ for method in {bcds,cxds,hybrid,scrublet}; do
    for data in {chcl,chpb,demu,hgmm}; do
      R --vanilla -f ./R/apply_method.R  --args --data=$data --method=$method
    done
  done
```
will run four methods on all data sets.

Results are written out in text files named according to ```./results/tmp/[dataset]_[method].txt```, which contain one line for each cell and two columns. The first
column contains the cell's barcode, the second the doublet score of the method.
Higher scores indicate more confidence in annotating a doublet.
For ```DoubletDecon``` we annotate doublet calls instead of scores.

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

* ```bcds_score```
* ```cxds_score```
* ```dblCells_score```
* ```dblDecon_call```
* ```dblDetection_score```
* ```dblFinder_score```
* ```hybrid_score```
* ```libsize_score```
* ```numfeat_score```
* ```scrublet_score```

This concludes doublet annotation. Processed and doublet-annotated data  can also be downloaded directly:

```{bash}
$ lftp ...
```

### Step 4: Analyzing doublet annotation results
We analyze computational doublet annoations in four ```SingleCellExperiment``` objects, one for each data set (see Step 3 above).

To reproduce tables and figures you can run the following scripts:

* ```./R/mk_tabs.R```
* ```./R/mk_figs.R```
* ```./R/mk_running-time.R```

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

#### Generate figures:

Figures are generated by the ```mk_figs.R``` script, analogously to the tables above:

```{bashr}
$ R --vanilla < ./R/mk_figs.R
```

 The figures generated are also written to the ```./results/``` directory. The following list of figures is produced:

* ```fig_cxds_[method].jpg```
* ```fig_perf_strat.pdf``` Methods performance for each data set.
* ```fig_comp.pdf```
* ```fig_size-strat_[method].pdf```
* ```fig_tpfn_[method].jpg```

#### Assess running time:

To assess running time, we limited container resources to two cores by using the following command:

```{bash}
$ docker run -it --cpuset-cpus="0-1"                                  \
                 --name scds_tab_fig                                  \
                 --net=host                                           \
                 --env DISPLAY=$DISPLAY                               \
                 --volume $HOME/.Xauthority:/home/rstudio/.Xauthority \
                 kostkalab_scds                                       \
                 /bin/bash
```

And then, after getting the code and data, unpacking and processing:

```{bash}
$ cd /home/rstudio/scds_manuscript
$ R --vanilla < ./R/mk_running-times.R
$ cat ./results/tab_running-times.tex
```
