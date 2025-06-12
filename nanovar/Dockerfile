FROM continuumio/miniconda3

ENV DEBIAN_FRONTEND=noninteractive

# Create group and user with UID/GID 1000
RUN groupadd -g 1000 nanovargroup && \
    useradd -m -u 1000 -g nanovargroup -s /bin/bash nanovar

# Set working directory
WORKDIR /app/nanovar

# Install conda environment with system bio tools
RUN conda create -n nanovar-env -c conda-forge -y \
    python=3.11 \
    bedtools=2.26.0 \
    samtools=1.3.0 \
    minimap2=2.17 && \
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
