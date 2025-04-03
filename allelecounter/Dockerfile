# Use a base image with Mamba
FROM mambaorg/micromamba:latest

# Set working directory
WORKDIR /opt/env

# Define the conda environment directly in the Dockerfile
RUN echo "name: allelecount_env\n\
channels:\n\
  - conda-forge\n\
  - bioconda\n\
  - defaults\n\
dependencies:\n\
  - allelecounter\n\
  - mamba" > /tmp/conda.yml

# Use Mamba to create the environment
RUN micromamba create -f /tmp/conda.yml -y && micromamba clean --all --yes

# Set the environment path to the created environment
ENV PATH="/opt/conda/envs/allelecount_env/bin:$PATH"

# Verify installation
RUN alleleCounter --version

# Set the default working directory for the container
WORKDIR /data

# Set the default command for the container
CMD ["/bin/bash"]
