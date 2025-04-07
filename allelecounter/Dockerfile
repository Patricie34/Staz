FROM continuumio/miniconda3

# Set environment variables
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

# Install cancerit-allelecount and its dependencies using conda
RUN conda install -c bioconda -c conda-forge -y \
    cancerit-allelecount \
    htslib \
    procps \
    && conda clean -afy \
    && find /opt/conda/ -follow -type f -name '*.a' -delete \
    && find /opt/conda/ -follow -type f -name '*.pyc' -delete \
    && find /opt/conda/ -follow -type f -name '*.js.map' -delete

# Set the working directory in the container
WORKDIR /app

# Command to run when the container starts
CMD ["/bin/bash"]
