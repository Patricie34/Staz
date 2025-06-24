# Base image with Miniconda3
FROM continuumio/miniconda3

# Set working directory
WORKDIR /app

# Copy environment.yml into the container
COPY environment.yml .

# Create the Conda environment
RUN conda env create -f environment.yml && \
    conda clean -a

# Activate the environment by adding it to PATH
ENV PATH /opt/conda/envs/savana/bin:$PATH

# Install Savana from GitHub, latest version
RUN git clone https://github.com/cortes-ciriano-lab/savana.git && \
    cd savana && \
    pip install . && \
    cd .. && rm -rf savana

# Create a non-root user
RUN useradd -m -u 1000 -s /bin/bash user && \
    chown -R user:user /app

# Switch to non-root user
USER 1000

# Default command
CMD ["bash"]
