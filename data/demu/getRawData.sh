#!/bin/bash

set -e
set -u

cd ./data/demu/
mkdir ./proc
mkdir ./raw 

#- NOTE: we only get the "2.1" cells. These should be the "control" cells without integrin stimulation:
#        These are downloaded in the processRawData.R script. 

#- Here we get singlet/doublet calls from the demuxlet github repository: 
wget --no-verbose --timestamping -o timestamp_calls.txt https://raw.githubusercontent.com/yelabucsf/demuxlet_paper_code/master/fig3/ye1.ctrl.8.10.sm.best

mv ye1.ctrl.8.10.sm.best ./raw/

#- check md5 sums
md5sum ./raw/* >./downloaded_md5sums.txt

cd ../..


