# Use the official Ubuntu image as the base
FROM ubuntu:20.04

# Install system dependencies
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
    software-properties-common \
    && rm -rf /var/lib/apt/lists/*

# Add CRAN repository and install R base version 4.4.0
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9 \
    && add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/' \
    && apt-get update \
    && apt-get install -y --no-install-recommends r-base=4.4.0-1.2004.0 \
    && rm -rf /var/lib/apt/lists/*

# Install R packages from CRAN
RUN R -e "install.packages(c('ggplot2', 'scales', 'gtable', 'grid'), repos='http://cran.rstudio.com/')"

# Install Bioconductor and DNAcopy package
RUN R -e "install.packages('BiocManager', repos='http://cran.rstudio.com/')" \
    && R -e "BiocManager::install('DNAcopy', update = TRUE, ask = FALSE)"

# Set the working directory (optional)
WORKDIR /app
