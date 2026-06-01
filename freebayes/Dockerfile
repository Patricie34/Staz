# =============================================================================
# Dockerfile for freebayes-filter
# Script to filter variant files produced by the FreeBayes variant caller
# Code is mostly copied from the bcbio package
# =============================================================================

FROM continuumio/miniconda3:latest

# Metadata
LABEL maintainer="medgen-bioinf"
LABEL description="Docker image for freebayes-filter with bcftools 1.9"
LABEL repository="https://gitlab.ics.muni.cz/medgen-bioinf/freebayes-filter"

# Prevent conda auto-updates and ensure non-interactive behavior
ENV CONDA_AUTO_UPDATE_CONDA=false
ENV CONDA_ALWAYS_YES=true

# =============================================================================
# Step 1: Install system dependencies required for bioinformatics tools
# =============================================================================
RUN apt-get update && apt-get install -y \
    gcc \
    g++ \
    make \
    cmake \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-gnutls-dev \
    libssl-dev \
    libncurses5-dev \
    git \
    wget \
    && rm -rf /var/lib/apt/lists/*

# =============================================================================
# Step 2: Configure conda channels (order matters for package resolution)
# =============================================================================
RUN conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda config --set channel_priority strict

# =============================================================================
# Step 3: Create conda environment with ALL bio-dependencies first
# Note: bcftools=1.9 is pinned as per Oct 2023 recommendations
# =============================================================================
RUN conda create -n freebayes_filter -c bioconda \
    samtools \
    vcflib \
    vt \
    bcftools=1.9 \
    tabix \
    python=3.9 \
    pip \
    -y

# =============================================================================
# Step 4: Create non-root user (UID/GID 1000 for typical host compatibility)
# =============================================================================
RUN groupadd -g 1000 user && \
    useradd -m -u 1000 -g 1000 -s /bin/bash user

# =============================================================================
# Step 5: Clone the repository and install the Python package via pip
# IMPORTANT: Conda deps first, then pip (as per best practices)
# =============================================================================
WORKDIR /app

RUN git clone https://gitlab.ics.muni.cz/medgen-bioinf/freebayes-filter.git .

# Install build tools and the bcbioloot_package
# Using conda run to ensure we're in the correct environment
RUN conda run -n freebayes_filter pip install --upgrade build && \
    conda run -n freebayes_filter python -m build bcbioloot_package/ && \
    conda run -n freebayes_filter pip install bcbioloot_package/

# =============================================================================
# Step 6: Set ownership and permissions for the user
# =============================================================================
RUN chown -R user:user /app && \
    chown -R user:user /opt/conda/envs/freebayes_filter

# =============================================================================
# Step 7: Set up working directory and entrypoint
# =============================================================================
WORKDIR /data

# Make the script accessible from anywhere
RUN ln -s /app/freebayes.py /opt/conda/envs/freebayes_filter/bin/freebayes-filter && \
    chown -R user:user /data

# Switch to non-root user
USER user

# Entry point allows running: docker run <image> -i input.vcf -r ref.fasta
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "freebayes_filter", "python3", "/app/freebayes.py"]
CMD ["--help"]
