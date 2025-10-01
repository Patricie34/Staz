FROM ubuntu:22.04

LABEL description="Simple Docker image with minimap2 and samtools"

# Install dependencies
RUN apt-get update && apt-get install -y \
    git \
    curl \
    make \
    gcc \
    g++ \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    wget \
    tar \
    bzip2 \
    xz-utils \
    && rm -rf /var/lib/apt/lists/*

# Install Samtools (includes htslib)
RUN curl -L -o samtools-1.20.tar.bz2 https://github.com/samtools/samtools/releases/download/1.20/samtools-1.20.tar.bz2 \
    && tar -xjf samtools-1.20.tar.bz2 \
    && cd samtools-1.20 \
    && ./configure --prefix=/usr/local \
    && make -j"$(nproc)" && make install \
    && cd .. && rm -rf samtools-1.20 samtools-1.20.tar.bz2

# Install Minimap2 v2.30
RUN curl -L -o minimap2-2.30.tar.bz2 https://github.com/lh3/minimap2/releases/download/v2.30/minimap2-2.30.tar.bz2 \
    && tar -xjf minimap2-2.30.tar.bz2 \
    && cd minimap2-2.30 \
    && make -j"$(nproc)" \
    && cp minimap2 /usr/local/bin/ \
    && cd .. && rm -rf minimap2-2.30 minimap2-2.30.tar.bz2

# Create non-root user
RUN useradd -m -u 1000 -s /bin/bash user

# Set working directory
WORKDIR /app
RUN chown user:user /app

# Switch to non-root user
USER 1000

ENTRYPOINT ["/bin/bash"]
