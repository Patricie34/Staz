# Use the official R base image
FROM r-base:4.0.5

# Install necessary system dependencies
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libxml2-dev \
    libssl-dev \
    libxt-dev

# Install required R packages from CRAN
RUN R -e "install.packages(c('ggplot2', 'scales', 'gtable', 'grid'), repos='http://cran.rstudio.com/')"

# Install Bioconductor and DNAcopy package
RUN R -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager'); BiocManager::install('DNAcopy')"

# Set the working directory
WORKDIR /app

