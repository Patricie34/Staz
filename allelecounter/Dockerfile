FROM debian:latest

# Install system utilities and dependencies
RUN apt-get update && apt-get install -y \
    curl ca-certificates bzip2 procps && \
    rm -rf /var/lib/apt/lists/*

# Install Micromamba
RUN curl -L https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj -C /usr/local

# Set up environment variables
ENV PATH=/usr/local/bin/micromamba/bin:$PATH
ENV MAMBA_ROOT_PREFIX=/opt/micromamba

# Create and activate a Micromamba environment
COPY conda.yml /tmp/conda.yml
RUN micromamba create -y -n allelecount_env -f /tmp/conda.yml

# Activate the environment
ENV PATH=/opt/micromamba/envs/allelecount_env/bin:$PATH

# Set working directory
WORKDIR /app

CMD ["/bin/bash"]
