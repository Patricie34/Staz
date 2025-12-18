FROM python:3.10-slim

# Install system dependencies
RUN apt-get update && apt-get install -y \
    gcc \
    g++ \
    && rm -rf /var/lib/apt/lists/*

# Install Python packages
RUN pip install --no-cache-dir \
    pandas \
    matplotlib \
    seaborn \
    scipy \
    numpy

# Create non-root user
RUN useradd -m -u 1000 -s /bin/bash user

# Set a working directory
WORKDIR /data
RUN chown -R user:user /data

# Switch to non-root user
USER 1000

# Default command
CMD ["python3"]
