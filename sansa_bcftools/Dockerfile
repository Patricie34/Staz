FROM continuumio/miniconda3

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    g++ \
    zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# Create group and user with UID 1000
RUN groupadd sansagroup && \
    useradd -m -u 1000 -g sansagroup -s /bin/bash sansa

# Set working directory
WORKDIR /app/sansa

# Install conda environment with system bio tools
RUN conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda create -n sansa-env -y \
        sansa=0.2.5 \
        htslib \
        bcftools \
&& conda clean -afy


# Update PATH to use the conda environment
ENV PATH=/opt/conda/envs/sansa-env/bin:$PATH

# Set permissions
RUN chown -R sansa:sansagroup /app/sansa

# Switch to non-root user
USER 1000

# Confirm working directory
WORKDIR /app/sansa
