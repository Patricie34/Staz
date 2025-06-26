# Base image with correct Python version
FROM python:3.10.15-slim

# Avoid interactive prompts
ENV DEBIAN_FRONTEND=noninteractive

# Install required system packages
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    curl \
    zlib1g-dev \
    procps \  
    gcc \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Install Sniffles2 and updated Python dependencies
RUN pip install --no-cache-dir --upgrade pip \
    && pip install sniffles==2.6.2 \
                 pysam==0.23.3 \
                 edlib>=1.3.9 \
                 psutil>=5.9.4

# Create non-root user with UID 1000
RUN useradd -m -u 1000 -s /bin/bash user

# Set working directory and fix permissions
WORKDIR /app
RUN chown user:user /app

# Switch to the created user
USER 1000

