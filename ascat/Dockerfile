# Use an official Rocker image for R
FROM rocker/tidyverse:4.2.3

# Install system dependencies for ASCAT and Bioconductor
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    gfortran \
    libgit2-dev \
    python3 \
    python3-pip \
    && rm -rf /var/lib/apt/lists/*

# Install Bioconductor packages (GenomicRanges & IRanges)
RUN R -e "install.packages('BiocManager')" && \
    R -e "BiocManager::install(c('GenomicRanges', 'IRanges'))"

# Install devtools to install ASCAT from GitHub
RUN R -e "install.packages('devtools')" && \
    R -e "devtools::install_github('VanLoo-lab/ascat/ASCAT')"

# Set the working directory in the container
WORKDIR /app
