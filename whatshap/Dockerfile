# Use official Ubuntu base image
FROM ubuntu:22.04

# Avoid interactive prompts during build
ENV DEBIAN_FRONTEND=noninteractive

# Install required system packages
RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 \
    python3-pip \
    python3-dev \
    build-essential \
    zlib1g-dev \
    htslib-tools \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Install WhatsHap using pip
RUN pip3 install --no-cache-dir --upgrade pip whatshap

# Create a non-root user with UID 1000 (group will be auto-created)
RUN useradd -m -u 1000 -s /bin/bash user

# Set working directory
WORKDIR /app

# Ensure the directory is owned by the new user
RUN chown user:user /app

# Switch to the created non-root user
USER 1000

# Default command
CMD ["bash"]
