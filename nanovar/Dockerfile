FROM continuumio/miniconda3

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    g++ \
    && rm -rf /var/lib/apt/lists/*

# Create group and user with UID 1000
RUN groupadd nanovargroup && \
    useradd -m -u 1000 -g nanovargroup -s /bin/bash nanovar

# Set working directory
WORKDIR /app/nanovar

# Install conda environment with system bio tools
RUN conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda create -n nanovar-env -y \
        python=3.11 \
        bedtools=2.31.1 \
        samtools=1.22 \
        minimap2=2.29 && \
    conda clean -afy

# Copy Python requirements
COPY requirements.txt .

# Install Python dependencies inside the conda environment
RUN /bin/bash -c "source activate nanovar-env && pip install --no-cache-dir -r requirements.txt"

# Update PATH to use the conda environment
ENV PATH /opt/conda/envs/nanovar-env/bin:$PATH

# Set permissions
RUN chown -R nanovar:nanovargroup /app/nanovar

# Switch to non-root user
USER 1000

# Confirm working directory
WORKDIR /app/nanovar
