FROM continuumio/miniconda3

# Install system dependencies
RUN apt-get update && apt-get install -y \
    procps \
    libhts-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install cancerit-allelecount and its dependencies
RUN conda install -c bioconda -c conda-forge -y \
    cancerit-allelecount \
    htslib \
    && conda clean -afy \
    && find /opt/conda/ -follow -type f -name '*.a' -delete \
    && find /opt/conda/ -follow -type f -name '*.pyc' -delete \
    && find /opt/conda/ -follow -type f -name '*.js.map' -delete

# Set the working directory in the container
WORKDIR /app

# Command to run when the container starts
CMD ["/bin/bash"]
