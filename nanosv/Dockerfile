FROM python:3.10.15-slim

ENV DEBIAN_FRONTEND=noninteractive

# Install NanoSV from PyPI
RUN pip install NanoSV

# Create non-root user with UID 1000
RUN useradd -m -u 1000 -s /bin/bash user

# Set working directory and fix permissions
WORKDIR /app
RUN chown user:user /app

# Switch to the created user
USER 1000
