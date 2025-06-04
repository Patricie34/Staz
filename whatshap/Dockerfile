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
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Install WhatsHap using pip
RUN pip3 install --no-cache-dir --upgrade pip whatshap 

# Create a group with GID 1 and user with UID 1000
RUN groupadd -g 1 group && \
    useradd -m -u 1000 -g 1 -s /bin/bash user

# Set working directory
WORKDIR /app

# Ensure the directory is owned by the new user
RUN chown user:group /app

# Switch to the created non-root user (UID 1000, GID 1)
USER 1000:1

# Default command
CMD ["bash"]
