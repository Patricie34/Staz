# Use an official Rocker image for R
FROM rocker/tidyverse:4.2.3 as builder

# Install system dependencies for ASCAT, Bioconductor, and alleleCount compilation
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    gfortran \
    libgit2-dev \
    python3 \
    python3-pip \
    build-essential \
    git \
    libboost-dev \
    libboost-iostreams-dev \
    libboost-filesystem-dev \
    && rm -rf /var/lib/apt/lists/*

# Install Bioconductor packages (GenomicRanges & IRanges)
RUN R -e "install.packages('BiocManager')" && \
    R -e "BiocManager::install(c('GenomicRanges', 'IRanges', 'DNAcopy'))"

# Install devtools to install ASCAT from GitHub
RUN R -e "install.packages('devtools')" && \
    R -e "devtools::install_github('VanLoo-lab/ascat/ASCAT')"

# Clone the alleleCount repository and build it
RUN git clone https://github.com/cancerit/alleleCount.git /alleleCount
WORKDIR /alleleCount

# Ensure all necessary build dependencies are installed before compiling
RUN apt-get update && apt-get install -y \
    libgsl-dev \
    libbz2-dev \
    liblzma-dev

# Try building the alleleCount tool
RUN make

# Copy the compiled alleleCount to the final image
FROM ubuntu:20.04

# Install required dependencies for alleleCount and ASCAT
RUN apt-get update && apt-get install -y \
    apt-transport-https \
    locales \
    curl \
    ca-certificates \
    libperlio-gzip-perl \
    bzip2 \
    psmisc \
    time \
    zlib1g \
    liblzma5 \
    libncurses5 \
    p11-kit \
    python3 \
    python3-pip \
    && rm -rf /var/lib/apt/lists/*

# Set environment variables for ASCAT and alleleCount
ENV OPT /opt/wtsi-cgp
ENV PATH=$OPT/bin:$OPT/biobambam2/bin:$PATH
ENV PERL5LIB=$OPT/lib/perl5
ENV LD_LIBRARY_PATH=$OPT/lib
ENV LC_ALL=en_US.UTF-8
ENV LANG=en_US.UTF-8

# Copy the compiled alleleCount from the builder stage
COPY --from=builder /alleleCount /opt/wtsi-cgp

# Set the working directory in the container
WORKDIR /app
