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
        autoconf \
        automake \
        libtool \
    && rm -rf /var/lib/apt/lists/*

# Install Bioconductor packages (GenomicRanges, IRanges, DNAcopy)
RUN R -e "install.packages('BiocManager')" && \
    R -e "BiocManager::install(c('GenomicRanges', 'IRanges', 'DNAcopy'))"

# Install devtools to install ASCAT from GitHub
RUN R -e "install.packages('devtools')" && \
    R -e "devtools::install_github('VanLoo-lab/ascat/ASCAT')"

# Clone the alleleCount repository from the 'dev' branch and build it
RUN git clone --branch dev https://github.com/cancerit/alleleCount.git /alleleCount && \
    cd /alleleCount/build && \
    ./build.sh

# Set environment variables for alleleCount (if needed)
ENV ALLELECOUNTPATH=/alleleCount/c

# Optionally, you can add your working directory here
WORKDIR /app

# Set the default command to start an R session (or run a script)
CMD ["R"]
