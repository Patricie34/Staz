FROM python:3.11-slim

# Install bedtools and standard utilities
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        bedtools \
        procps \
        sed \
        coreutils \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Create non-root user
RUN useradd -m -u 1000 -s /bin/bash user

# Set a working directory
WORKDIR /data
RUN chown -R user:user /data

# Switch to non-root user
USER 1000

# Default command
CMD ["bash"]
