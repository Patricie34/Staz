# Use a base image with Mamba
FROM mambaorg/micromamba:latest

# Set working directory for environment creation
WORKDIR /opt/env

# Define the conda environment file in the Dockerfile
COPY conda.yml /tmp/conda.yml

# Use Mamba to create the environment
RUN micromamba create -f /tmp/conda.yml -y && micromamba clean --all --yes

# Activate the environment and set the path
ENV PATH="/opt/conda/envs/allelecount_env/bin:$PATH"

# Verify the installation of alleleCounter
RUN alleleCounter --version

# Set the default working directory for the container
WORKDIR /app

# Set the default command for the container to launch a shell
CMD ["/bin/bash"]
