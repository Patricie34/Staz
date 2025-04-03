# Use the Micromamba base image
FROM mambaorg/micromamba:latest

# Install required system utilities (including 'ps' for Nextflow)
RUN apt-get update && apt-get install -y procps

# Set up your Micromamba environment
COPY environment.yml /tmp/environment.yml
RUN micromamba create -y -n myenv -f /tmp/environment.yml

# Activate the environment when running the container
ENV MAMBA_DEFAULT_ENV=myenv
ENV PATH=/opt/conda/envs/myenv/bin:$PATH

# Set the working directory
WORKDIR /app

# Define the entrypoint (optional)
CMD ["/bin/bash"]
