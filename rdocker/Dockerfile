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
WORKDIR /usr/local/src/myscripts

# Copy your script to the Docker image (replace your_script.R with your actual script name)
COPY your_script.R .

# Default command to run your script (replace your_script.R with your actual script name)
CMD ["Rscript", "your_script.R"]
