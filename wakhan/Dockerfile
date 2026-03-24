FROM continuumio/miniconda3

# Set working directory
WORKDIR /app

# Install system dependencies (root by default - no USER root needed)
RUN apt-get update && apt-get install -y \
    git \
    && rm -rf /var/lib/apt/lists/*

# Copy environment first (better layer caching)
COPY environment.yml .

# Create the conda environment
RUN conda env create -f environment.yml -n wakhan && \
    conda clean -a

# Add conda env to PATH
ENV PATH /opt/conda/envs/wakhan/bin:$PATH

# Clone Wakhan from source
RUN git clone https://github.com/KolmogorovLab/Wakhan.git /app/Wakhan

WORKDIR /app/Wakhan

# Create a non-root user with UID 1000
RUN useradd -m -u 1000 -s /bin/bash user && \
    chown -R user:user /app

# Switch to non-root user
USER 1000

# Install Wakhan package
RUN python -m pip install --no-cache-dir -e .

# Default command
CMD ["wakhan", "--help"]
