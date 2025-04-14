# Use base image with Miniconda
FROM continuumio/miniconda3:latest

# Install system packages needed for devtools and R compilation
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

# Copy Conda environment file
COPY ascat_env.yml /app/ascat_env.yml

# Create Conda environment
RUN conda env create -f /app/ascat_env.yml --prefix /opt/conda/envs/ascat_env

# Initialize Conda for bash shell
RUN /opt/conda/bin/conda init bash

# Set the environment variables to ensure Conda and R are available
ENV PATH="/opt/conda/envs/ascat_env/bin:$PATH"
ENV CONDA_PREFIX="/opt/conda/envs/ascat_env"
ENV CONDA_DEFAULT_ENV="ascat_env"

# Install ASCAT from GitHub using Rscript from the environment
RUN /opt/conda/envs/ascat_env/bin/Rscript -e "devtools::install_github('VanLoo-lab/ascat/ASCAT')"

# Set default command
CMD ["bash"]
