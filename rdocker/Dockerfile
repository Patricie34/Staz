# Use the official Ubuntu image as the base
FROM ubuntu:20.04

# Install system dependencies and R
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libxml2-dev \
    libssl-dev \
    libxt-dev \
    liblzma-dev \
    zlib1g-dev \
    libgmp-dev \
    libgit2-dev \
    gfortran \
    && apt-get install -y --no-install-recommends r-base \
    && rm -rf /var/lib/apt/lists/*

# Install R packages from CRAN and Bioconductor
RUN R -e "install.packages(c('ggplot2', 'scales', 'gtable', 'grid'), repos='http://cran.rstudio.com/')" \
    && R -e "if (!require('BiocManager', quietly = TRUE)) install.packages('BiocManager'); BiocManager::install('DNAcopy', update = TRUE, ask = FALSE)"

# Set the working directory (optional)
WORKDIR /app
