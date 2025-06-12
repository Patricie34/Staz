# Use an official minimal Python runtime as a parent image
FROM python:3.9-slim

# Update package lists and install necessary system dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends git build-essential && \
    rm -rf /var/lib/apt/lists/*

# Create non-root user with UID 1000 and add to the group with GID 1
RUN addgroup --gid 1 cutesvgroup && \
    useradd -m -u 1000 -s /bin/bash -g cutesvgroup cutesvuser

# Set working directory and fix permissions
WORKDIR /app
RUN chown -R cutesvuser:cutesvgroup /app

# Switch to the created user
USER 1000

# Clone the cuteSV repository
RUN git clone https://github.com/tjiangHIT/cuteSV.git /app/cuteSV

# Install cuteSV using setup.py
WORKDIR /app/cuteSV
RUN python setup.py install

