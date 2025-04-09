# Use a base image with Miniconda3
FROM continuumio/miniconda3:latest

# Set the working directory
WORKDIR /app

# Install dependencies from conda YML conda
COPY conda.yml /app/conda.yml

# Install the conda via conda
RUN conda env create -f conda.yml

# Set the default conda to activate
SHELL ["conda", "run", "-n", "base", "/bin/bash", "-c"]

# Set the default command for the container
CMD ["bash"]
