FROM python:3.8-slim

ENV DEBIAN_FRONTEND=noninteractive

# Install NanoSV from PyPI
RUN pip install NanoSV

# Create non-root user
RUN useradd -m -u 1000 -s /bin/bash user

# Set working directory and fix permissions
WORKDIR /app
RUN chown user:user /app

# Switch to the created user
USER 1000
