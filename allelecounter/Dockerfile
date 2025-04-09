# Use a base image with Miniconda3
FROM continuumio/miniconda3:latest

# Set the working directory
WORKDIR /app

# Copy the environment YAML file
COPY conda.yml /app/conda.yml

# Create the environment with a custom path
RUN conda env create -f /app/conda.yml --prefix /opt/conda/envs/myenv

# Activate the environment
RUN echo "conda activate myenv" >> ~/.bashrc

# Set the default environment to activate
SHELL ["conda", "run", "-n", "myenv", "/bin/bash", "-c"]

# Set the default command for the container
CMD ["bash"]
