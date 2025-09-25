# Use the existing Pod5
FROM quay.io/biocontainers/pod5:0.3.33--pyhdfd78af_0

# Install rsync
RUN apt-get update && \
    apt-get install -y rsync && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Create a non-root user with UID 1000 
RUN useradd -m -u 1000 -s /bin/bash user

# Set working directory
WORKDIR /app

# Ensure directory is owned by the new user
RUN chown user:user /app

# Switch to non-root user
USER 1000
