# Use official Ubuntu base image
FROM ubuntu:24.04

# Avoid interactive prompts during build
ENV DEBIAN_FRONTEND=noninteractive

# Install required system packages
RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 \
    python3-pip \
    python3-dev \
    build-essential \
    zlib1g-dev \
    curl \
    make \
    gcc \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Install WhatsHap using pip
RUN pip3 install --no-cache-dir --upgrade pip whatshap

# Install HTSlib (to get bgzip and tabix)
ENV HTSLIB_VERSION=1.19
RUN curl -fsSL https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 | tar -xj && \
    cd htslib-${HTSLIB_VERSION} && \
    ./configure --prefix=/usr/local && \
    make -j$(nproc) && \
    make install && \
    cd .. && rm -rf htslib-${HTSLIB_VERSION

# Create a non-root user with UID 1000 (group will be auto-created)
RUN useradd -m -u 1000 -s /bin/bash user

# Set working directory
WORKDIR /app

# Ensure the directory is owned by the new user
RUN chown user:user /app

# Switch to the created non-root user
USER 1000

# Default command
CMD ["bash"]
