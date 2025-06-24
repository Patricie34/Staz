# Base image with Miniconda and Python 3.8+
FROM continuumio/miniconda3

# Set working directory for environment setup
WORKDIR /app

# Copy environment configuration
COPY environment.yml .

# Create the Conda environment
RUN conda env create -f environment.yml && \
    conda clean -a

# Add Conda env to PATH
ENV PATH /opt/conda/envs/severus/bin:$PATH

# Install Severus from GitHub
RUN git clone https://github.com/KolmogorovLab/Severus.git && \
    cd Severus && \
    python setup.py install && \
    cd .. && rm -rf Severus

# Create a non-root user with UID 1000
RUN useradd -m -u 1000 -s /bin/bash user && \
    chown -R user:user /app

# Switch to non-root user
USER 1000

# Default shell
CMD ["bash"]
