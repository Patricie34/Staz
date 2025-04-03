# Use a base image with Mamba
FROM mambaorg/micromamba:latest

# Set working directory
WORKDIR /opt/env

# Copy conda.yml into the container
COPY conda.yml /tmp/conda.yml

# Use Mamba to create the environment
RUN micromamba create -f /tmp/conda.yml -y && micromamba clean --all --yes

# Set the environment path to the created environment
ENV PATH="/opt/conda/envs/allelecount_env/bin:$PATH"

# Set the default working directory for the container
WORKDIR /app

# Set the default command for the container
CMD ["/bin/bash"]
