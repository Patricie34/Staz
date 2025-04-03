FROM mambaorg/micromamba:1.5.4-debian-bullseye

# Install system utilities (procps includes 'ps' for Nextflow)
RUN apt-get update && apt-get install -y procps

# Set up your Micromamba environment
COPY conda.yml /tmp/conda.yml
RUN micromamba create -y -n allelecounter-env -f /tmp/conda.yml

# Activate the environment when running the container
ENV MAMBA_DEFAULT_ENV=allelecounter-env
ENV PATH=/opt/conda/envs/allelecounter-env/bin:$PATH

# Set the working directory
WORKDIR /app

# Define the entrypoint (optional)
CMD ["/bin/bash"]
