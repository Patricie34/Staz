# Use an official Miniconda image as the base
FROM continuumio/miniconda3:latest

# Set environment variables
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

# Install essential utilities and cancerit-allelecount
RUN apt-get update && apt-get install -y procps && \
    conda install -c bioconda -c conda-forge -y cancerit-allelecount && \
    conda clean -afy && \
    find /opt/conda/ -follow -type f -name '*.a' -delete && \
    find /opt/conda/ -follow -type f -name '*.pyc' -delete && \
    find /opt/conda/ -follow -type f -name '*.js.map' -delete && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Set the working directory in the container
WORKDIR /app

# Command to run when the container starts
CMD ["/bin/bash"]
