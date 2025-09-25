# Start from a minimal Ubuntu base
FROM ubuntu:22.04

# Avoid interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Update and install required packages: python3, pip, rsync, bash, curl, etc.
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        python3 \
        python3-pip \
        rsync \
        bash \
        curl \
        git \
        && apt-get clean && rm -rf /var/lib/apt/lists/*

# Upgrade pip
RUN python3 -m pip install --upgrade pip

# Install latest pod5 via pip
RUN pip install --no-cache-dir --upgrade pod5

# Optional: create a non-root user
RUN useradd -m -u 1000 -s /bin/bash user

# Set working directory
WORKDIR /app
RUN chown user:user /app

# Switch to non-root user
USER 1000

# Default shell
SHELL ["/bin/bash", "-c"]

# Optional: enable pod5 debug by default
ENV POD5_DEBUG=1
