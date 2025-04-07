FROM debian:latest

# Install system utilities and dependencies
RUN apt-get update && apt-get install -y \
    curl ca-certificates bzip2 procps bash tar && \
    rm -rf /var/lib/apt/lists/*

# Install Micromamba
RUN curl -L https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj -C /usr/local

# Environment variables
ENV MAMBA_ROOT_PREFIX=/opt/micromamba
ENV PATH=/opt/micromamba/envs/allelecount_env/bin:/usr/local/bin/micromamba/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

# Copy environment definition
COPY conda.yml /tmp/conda.yml

# Create micromamba environment
RUN micromamba create -y -n allelecount_env -f /tmp/conda.yml && \
    micromamba clean --all --yes

# Set working directory
WORKDIR /app

# Make sure entrypoint doesn't override bash
ENTRYPOINT []
CMD ["/bin/bash"]
