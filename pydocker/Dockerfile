# Use an official Python base image
FROM python:3.9-slim

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive

# Install samtools, bcftools, and any dependencies
RUN apt-get update && \
    apt-get install -y \
    samtools \
    bcftools \
    && rm -rf /var/lib/apt/lists/*

# Set up a working directory
WORKDIR /usr/src/app


