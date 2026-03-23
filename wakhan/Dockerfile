FROM mambaorg/micromamba:1.5.8

# Install system dependencies
RUN apt-get update && apt-get install -y \
    git \
    && rm -rf /var/lib/apt/lists/*

# Copy environment first (better layer caching)
COPY environment.yml /tmp/environment.yml

# Create environment from source
RUN micromamba create -y -n wakhan -f /tmp/environment.yml && \
    micromamba clean --all --yes

# Set environment variables
ENV MAMBA_DEFAULT_ENV=wakhan \
    MAMBA_ROOT_PREFIX=/opt/conda \
    PATH=/opt/conda/envs/wakhan/bin:$PATH

# Clone Wakhan from source
RUN git clone https://github.com/KolmogorovLab/Wakhan.git /app/Wakhan

WORKDIR /app/Wakhan

# Create a non-root user with UID 1000
RUN useradd -m -u 1000 -s /bin/bash user && \
    chown -R user:user /app/Wakhan && \
    chown -R user:user /opt/conda/envs/wakhan

# Switch to non-root user
USER 1000

# Install Wakhan package
RUN python -m pip install --no-cache-dir -e .

# Default command
ENTRYPOINT ["python", "wakhan.py"]
CMD ["--help"]
