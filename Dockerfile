FROM r-base:4.4.1
ARG buildDate

# Labels
LABEL about.image="R base image to run CloEvoScore pipeline"
LABEL image.authors="Gaia Mazzocchetti, Andrea Poletti, Vincenza Solli"
LABEL image.maintainer="Gaia Mazzocchetti"
LABEL image.maintainer.email="gaia.mazzocchetti2@unibo.it"
LABEL base.image="r-base:4.4.1"
LABEL image.version="0.1"
LABEL buildDate=${buildDate}
LABEL software="CloEvoScore pipeline aim to estimate evolution patient trajectory"
LABEL software.version="1.0.0"
LABEL about.tags="Multiple Myeloma, Clonal Evolution, Copy Number Alteration"

# Install required system libraries
RUN apt-get update && apt-get install -y bash \
    libcurl4-openssl-dev \
    libssl-dev \
    libssh2-1-dev \
    procps \
    libxml2-dev \
    libpng-dev \
    libnlopt-dev \
    libfreetype6-dev \
    libfontconfig1-dev \
    libjpeg-dev && \
    apt-get clean && rm -rf /var/lib/apt/lists/*


# Install R dependencies in one step to optimize caching
RUN R -e "options(repos = c(CRAN = 'https://cran.r-project.org')); install.packages(c('BiocManager','factoextra','purrr', 'optparse', 'data.table', 'dplyr', 'stringr', 'ggplot2', 'readr', 'dbscan', 'proxy', 'ggrepel', 'ggforce'), dependencies=TRUE)"

# Install Bioconductor dependencies
RUN R -e "BiocManager::install(c('GenomicRanges', 'plyranges'))"

# Set working directory
WORKDIR /cloevoscore

# Copy scripts and data
COPY script /cloevoscore/script/
COPY www/data /cloevoscore/data/

# Ensure the main script is executable
RUN chmod +x /cloevoscore/script/Main_script.R

## Set entry point but allow arguments to be passed
ENTRYPOINT ["Rscript", "/cloevoscore/script/Main_script.R"]

