# ------------------------------------------------------------
#  Dockerfile: samtools + minimap2 (latest versions)
# ------------------------------------------------------------
#   Base image
FROM ubuntu:22.04

#   Metadata
LABEL description="Ubuntu 22.04 image with the latest Samtools and Minimap2 built from source."

# ------------------------------------------------------------
#   1. Install system‑level build dependencies
# ------------------------------------------------------------
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        ca-certificates \
        wget \
        curl \
        gnupg \
        build-essential \
        libncurses5-dev \
        libncursesw5-dev \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        libcurl4-gnutls-dev \
        libssl-dev \
        libgomp1 \
        python3 \
        python3-pip \
        git && \
    rm -rf /var/lib/apt/lists/*

# ------------------------------------------------------------
#   2. Install the latest Samtools (compiled from source)
# ------------------------------------------------------------
ARG SAMTOOLS_REPO="https://github.com/samtools/samtools.git"
ARG SAMTOOLS_TAG="master"          # change to a tag/branch if you want a specific release

RUN git clone --depth 1 --branch --recurse-submodules ${SAMTOOLS_TAG} ${SAMTOOLS_REPO} /opt/samtools && \
    cd /opt/samtools && \
    make -j$(nproc) && \
    make install && \
    cd / && \
    rm -rf /opt/samtools

# ------------------------------------------------------------
#   3. Install the latest Minimap2 (compiled from source)
# ------------------------------------------------------------
ARG MINIMAP2_REPO="https://github.com/lh3/minimap2.git"
ARG MINIMAP2_TAG="master"          # change to a tag/branch if you want a specific release

RUN git clone --depth 1 --branch ${MINIMAP2_TAG} ${MINIMAP2_REPO} /opt/minimap2 && \
    cd /opt/minimap2 && \
    make -j$(nproc) && \
    cp minimap2 mm2-fast mm2-precise /usr/local/bin/ && \
    cd / && \
    rm -rf /opt/minimap2

# ------------------------------------------------------------
#   4. Verify installations (optional but nice for debugging)
# ------------------------------------------------------------
RUN samtools --version && \
    minimap2 --version
    
# Optional: create a non-root user
RUN useradd -m -u 1000 -s /bin/bash user

# Set working directory
WORKDIR /app
RUN chown user:user /app

# Switch to non-root user
USER 1000
