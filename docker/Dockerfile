
FROM rocker/verse:3.5.2 

RUN apt-get update 							\
&& DEBIAN_FRONTEND=noninteractive apt-get -qq install libhdf5-dev	\
&& DEBIAN_FRONTEND=noninteractive apt-get -qq install bzip2		\
&& DEBIAN_FRONTEND=noninteractive apt-get -qq install libudunits2-dev

#- Install R packages 
#====================
RUN Rscript -e 'install.packages("BiocManager") ; BiocManager::install()' 		\
&&  Rscript -e 'library("BiocManager") ; BiocManager::install("GenomicRanges")'		\
&&  Rscript -e 'library("BiocManager") ; BiocManager::install("GenomicFeatures")'	\
&&  Rscript -e 'library("BiocManager") ; BiocManager::install("Seurat")'		\
&&  Rscript -e 'library("BiocManager") ; BiocManager::install("scran")'			\
&&  Rscript -e 'library("BiocManager") ; BiocManager::install("scater")'		\
&&  Rscript -e 'library("BiocManager") ; BiocManager::install("DropletUtils")'		\
&&  Rscript -e 'library("BiocManager") ; BiocManager::install("xgboost")'		\
&&  Rscript -e 'library("BiocManager") ; BiocManager::install("BiocFileCache")'		\
&&  Rscript -e 'library("BiocManager") ; BiocManager::install("PRROC")'			\
&&  Rscript -e 'library("BiocManager") ; BiocManager::install("pROC")'			\
&&  Rscript -e 'library("BiocManager") ; BiocManager::install("rsvd")'	

RUN Rscript -e 'devtools::install_github("chris-mcginnis-ucsf/DoubletFinder")'		\
&&  Rscript -e 'library("BiocManager") ; BiocManager::install("fields")'		\
&&  Rscript -e 'library("BiocManager") ; BiocManager::install("modes")'			\
&&  Rscript -e 'library("BiocManager") ; BiocManager::install("clusterProfiler")'	\
&&  Rscript -e 'library("BiocManager") ; BiocManager::install("hopach")'		\
&&  Rscript -e 'library("BiocManager") ; BiocManager::install("mygene")'		\
&&  Rscript -e 'library("BiocManager") ; BiocManager::install("MCL")'			\
&&  Rscript -e 'library("BiocManager") ; BiocManager::install("mygene")'		\
&&  Rscript -e 'library("BiocManager") ; BiocManager::install("DeconRNASeq")'		\
&&  Rscript -e 'library("BiocManager") ; BiocManager::install("kableExtra")'		\
&&  Rscript -e 'devtools::install_github("EDePasquale/DoubletDecon")'			\
&&  Rscript -e 'devtools::install_github("kostkalab/scds")'

#- Install miniconda for python 3
#================================
#- don't use $HOME, so here is where miniconda goes (MCL):
RUN mkdir /opt/user && chmod 755 /opt/user 				\
&&  mkdir /opt/user/rstudio && chown rstudio:rstudio /opt/user/rstudio  
ENV MCL="/opt/user/rstudio/miniconda"

USER rstudio

RUN cd ~										\
&& wget https://repo.continuum.io/miniconda/Miniconda3-3.19.0-Linux-x86_64.sh		\
&& bash Miniconda3-3.19.0-Linux-x86_64.sh -p ${MCL} -b				\
&& rm Miniconda3-3.19.0-Linux-x86_64.sh 
ENV PATH="${MCL}/bin:${PATH}"
RUN conda update -y conda								

#- Install scrublet and DoubletDetection
#=======================================
RUN  conda install -y python=3.7.2							\
&& conda update pip						

RUN cd ~										\
&& git clone https://github.com/AllonKleinLab/scrublet.git				\
&& cd scrublet										\
&& cat requirements.txt | sed '/annoy/d ; /umap-learn/d' > requirements_conda.txt	\
&& conda install -y --file ./requirements_conda.txt					\
&& yes w | pip install  annoy								\
&& yes w | pip install  umap-learn							\
&& yes w | pip install  --upgrade .							\
&& cd .. && rm -rf scrublet

RUN cd ~										 \
&& git clone https://github.com/JonathanShor/DoubletDetection.git			 \
&& cd DoubletDetection									 \
&& cat requirements.txt | sed '/git/d ; /tables/d ; /sklearn/d' > requirements_conda.txt \
&& conda install -y --file ./requirements_conda.txt					 \
&& conda install -y pytables 								 \
&& yes w | pip install git+https://github.com/DmitryUlyanov/Multicore-TSNE.git		 \
&& yes w | pip install git+https://github.com/JonathanShor/PhenoGraph.git		 \
&& yes w | pip install --upgrade .							 \
&&  cd .. && rm -rf DoubletDetection

ENV LD_PRELOAD="${MCL}/lib/libstdc++.so"
ENV LD_LIBRARY_PATH="${MCL}/lib:${LD_LIBRARY_PATH}"

RUN echo "source activate" > ~/.bashrc	



