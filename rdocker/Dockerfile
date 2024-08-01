# Use the rocker/r-ver base image for R version 4.4.0
FROM rocker/r-ver:4.4.0

# Install system dependencies for the packages
RUN apt-get update && \
    apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install the required CRAN packages
RUN R -e "install.packages(c('ggplot2', 'scales', 'gtable', 'grid'), repos='http://cran.rstudio.com/')"

# Install Bioconductor packages
RUN R -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager'); \
           BiocManager::install(version = '3.19'); \
           BiocManager::install('DNAcopy')"

# Set the working directory
WORKDIR /home/rstudio

