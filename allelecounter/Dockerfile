# Use a base image with Miniconda3
FROM continuumio/miniconda3:latest

# Install system dependencies needed for R packages and devtools
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    git \
    build-essential \
    procps \
    && apt-get clean

# Set working directory
WORKDIR /app

# Copy the conda environment file
COPY ascat_env.yml /app/ascat_env.yml

# Create the environment
RUN conda env create -f /app/ascat_env.yml --prefix /opt/conda/envs/ascat_env

# Set the environment to use in following RUN commands
SHELL ["conda", "run", "-n", "ascat_env", "/bin/bash", "-c"]

# Install ASCAT from GitHub (R package)
RUN Rscript -e "devtools::install_github('VanLoo-lab/ascat/ASCAT')"

# Default command
CMD ["bash"]
