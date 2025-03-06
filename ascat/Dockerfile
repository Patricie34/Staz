# 1️⃣ Stage 1: Compile alleleCounter
FROM ubuntu:20.04 AS builder

# Install dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    git \
    zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# Clone and compile alleleCounter
WORKDIR /opt
RUN git clone --recursive https://github.com/cancerit/alleleCount.git && \
    cd alleleCount && \
    make

# 2️⃣ Stage 2: Final Image with ASCAT + alleleCounter
FROM rocker/tidyverse:4.2.3

# Install additional system dependencies for ASCAT
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

# Copy compiled alleleCounter binary from builder
COPY --from=builder /opt/alleleCount/alleleCounter /usr/local/bin/

# Verify installation
RUN which alleleCounter && alleleCounter --version

# Set working directory
WORKDIR /app
