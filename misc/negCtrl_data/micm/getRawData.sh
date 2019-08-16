#!/bin/bash
############################################################
set -e
set -u

mkdir ./proc
mkdir ./raw 

############################################################
## Mohammed2017 Cell Reports https://www.sciencedirect.com/science/article/pii/S2211124717309610?via%3Dihub#app3
## GSE100597
############################################################
wget --no-verbose --timestamping -o timestamp_expression.txt \
     ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE100nnn/GSE100597/suppl/GSE100597%5Fcount%5Ftable%5FQC%5Ffiltered%2Etxt%2Egz

mv GSE100597_count_table_QC_filtered.txt.gz  ./raw/

############################################################
############################################################
############################################################
