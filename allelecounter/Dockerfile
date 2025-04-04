# Use a modern Ubuntu base image
FROM ubuntu:22.04  

# Install system dependencies
RUN apt-get update && apt-get install -y \
    curl ca-certificates procps git build-essential cmake \
    && rm -rf /var/lib/apt/lists/*

# Clone and install alleleCounter
WORKDIR /opt
RUN git clone https://github.com/cancerit/alleleCount.git && \
    cd alleleCount && \
    mkdir build && cd build && \
    cmake .. && make -j4 && \
    cp alleleCounter /usr/local/bin/

# Set working directory
WORKDIR /app

# Default command
CMD ["/bin/bash"]
