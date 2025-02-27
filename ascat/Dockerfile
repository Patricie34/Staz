# Use rocker/tidyverse as the base image
FROM rocker/tidyverse:4.2.3

# Set environment variables to avoid interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive \
    LC_ALL=C.UTF-8 \
    LANG=C.UTF-8

# Update and install system dependencies
RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y --no-install-recommends \
        build-essential \
        libcurl4-gnutls-dev \
        libxml2-dev \
        libssl-dev \
        wget \
        git \
        vim \
        g++ \
        make \
    && rm -rf /var/lib/apt/lists/*

# Install Bioconductor packages (GenomicRanges, IRanges, DNAcopy)
RUN R -e "install.packages('BiocManager')" && \
    R -e "BiocManager::install(c('GenomicRanges', 'IRanges', 'DNAcopy'))"

# Install devtools to install ASCAT from GitHub
RUN R -e "install.packages('devtools')" && \
    R -e "devtools::install_github('VanLoo-lab/ascat/ASCAT')"

# Clone the alleleCount repository and build it
RUN git clone https://github.com/cancerit/alleleCount.git /alleleCount && \
    cd /alleleCount && \
    make

# Set environment variables for alleleCount (if needed)
ENV ALLELECOUNTPATH=/alleleCount

# Optionally, you can add your working directory here
WORKDIR /app
