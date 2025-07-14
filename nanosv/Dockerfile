FROM python:3.8-slim

ENV DEBIAN_FRONTEND=noninteractive

# Install procps and samtools
RUN apt-get update && apt-get install -y --no-install-recommends \
    procps \
    samtools \
    && rm -rf /var/lib/apt/lists/*

# Install NanoSV from PyPI
RUN pip install NanoSV

# Create non-root user
RUN useradd -m -u 1000 -s /bin/bash user
WORKDIR /app
RUN chown user:user /app
USER 1000
