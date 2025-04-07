# Use an official Miniconda image as the base
FROM continuumio/miniconda3:latest

# Set maintainer label
LABEL maintainer="your-email@example.com"

# Set environment variables
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

# Install cancerit-allelecount and any other necessary tools
RUN conda install -c bioconda -c conda-forge -y cancerit-allelecount && \
    conda clean -afy && \
    find /opt/conda/ -follow -type f -name '*.a' -delete && \
    find /opt/conda/ -follow -type f -name '*.pyc' -delete && \
    find /opt/conda/ -follow -type f -name '*.js.map' -delete

# Set the working directory in the container
WORKDIR /app

# Command to run when the container starts
CMD ["/bin/bash"]
