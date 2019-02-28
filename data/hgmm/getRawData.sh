#!/bin/bash

set -e
set -u

#- need to download data via the browser:
#
# https://support.10xgenomics.com/single-cell-gene-expression/datasets/2.1.0/hgmm_12k
#
# Files downloaded:
#
# Gene / cell matrix HDF5 (raw) (md5 from website: 916c07c7659b06fc673fc6aeedba5693)
#
# Clustering analysis (md5 from website: 2540bb5b42f649f6266b6ef95cce0ef2)

cd ./data/hgmm

#- check md5 sums
md5sum ./raw/* >./downloaded_md5sums.txt

#- unpack cell annotations
cd ./raw && tar xvfz ./hgmm_12k_analysis.tar.gz

cd ../../..
