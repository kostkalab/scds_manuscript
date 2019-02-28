#!/bin/bash

set -e
set -u

#- THIS FOLLOWS THE INSTRUCTIONS ON
# https://satijalab.org/seurat/hashing_vignette.html
#
# - There are two cell-hashing datasets 
#
# (1) 8-HTO from human PBMC cells
#
# (2) 12-HTO dataset from four human cell lines
#
# we download both from Dropbox:
#
# https://www.dropbox.com/sh/c5gcjm35nglmvcv/AABGz9VO6gX9bVr5R2qahTZha?dl=0
#
# -> do directory download into ./raw
# -> unzip HTODemuxFiles.zip
# -> keep only the pbmc related files (compare with chcl)
#    - 
#    - 


#- check md5 sums
md5sum ./data/chpb/raw/* >./data/chpb/downloaded_md5sums.txt



###########################################################
###########################################################
