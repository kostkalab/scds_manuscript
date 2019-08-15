
#- Start off witht the v1 version of the docker image
docker run -it --name scdsR2 --net=host --env=DISPLAY --volume $HOME/.Xauthority:/home/rstudio/.Xauthority kostkalab/scds:v1.0  /bin/bash


#===================
#-  INSTALL NEW SCDS
#===================

#- Need the revisions for scds - use version 1.1.2,
#  commit edf0d5fd1924492f0450ad9e63dc7ae081171295
#   need to work arond the R 3.6 requirement...
cd
Rscript -e 'BiocManager::install("BiocStyle")' #- required
git clone -n https://github.com/kostkalab/scds.git && cd ./scds
git checkout edf0d5fd1924492f0450ad9e63dc7ae081171295
cat DESCRIPTION  | sed 's/Depends: R (>= 3.6.0)/Depends: R (>= 3.5.2)/' > ./tmp && mv ./tmp DESCRIPTION
cd .. && R CMD build ./scds && R CMD INSTALL ./scds_1.1.2.tar.gz

#=====================
#- INSTALL additional packages
#=====================
Rscript -e 'BiocManager::install("sctransform")'
Rscript -e 'BiocManager::install("pamr")'

exit

#- RE-START
docker start scdsR2

#- make an actual image locally:
DI=$(docker ps | grep scdsR2 | cut -d " " -f 1)
docker commit $DI kostkalab/scds #- new  image with changes

#- log into docker hub
docker login --username=***** --email=*****
#- fish out the image  ID using docker images and looking at time created and repository
docker tag 07b14efc829a kostkalab/scds:v2.0
#- publish
docker push kostkalab/scds:v2.0
