FROM r-base:4.4.1
ARG buildDate

COPY script/ cloevoscore/

LABEL about.image="R base image to run CloEvoScore pipeline"
LABEL image.authors="Gaia Mazzocchetti, Andrea Poletti, Vincenza Solli"
LABEL image.maintainer="Gaia Mazzocchetti"
LABEL image.maintainer.email="gaia.mazzocchetti2@unibo.it"
LABEL base.image="r-base:4.4.1"
LABEL image.version="0.1"
LABEL buildDate=${buildDate}
LABEL software="CloEvoScore pipeline aim to estimate evolution patient trajectory"
LABEL software.version="0.2.0"
LABEL about.tags="Multiple Myeloma, Clonal Evolution"

RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libssh2-1-dev \
    procps \
    libxml2-dev \
    libpng-dev \
    libjpeg-dev
	 
RUN R -e "install.packages(c('BiocManager', 'optparse','data.table','dplyr','stringr','ggplot2','readr', 'dbscan', 'factoextra', 'proxy', 'ggrepel', 'ggforce' ), repos='http://cran.rstudio.com/')"
RUN R -e "BiocManager::install(c( 'GenomicRanges'))"

ENTRYPOINT ["Rscript"]
