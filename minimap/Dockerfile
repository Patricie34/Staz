# -------- Base image --------
FROM ubuntu:22.04

LABEL maintainer="your_email@example.com" \
      description="Docker image with minimap2 v2.30, samtools v1.20, htslib v1.20, non-root user"

# -------- Install build dependencies --------
RUN apt-get update && \
    apt-get install -y \
        git \
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
        pkg-config \
        && rm -rf /var/lib/apt/lists/*

# -------- Set versions and build directory --------
ENV HTSLIB_VERSION=1.20 \
    SAMTOOLS_VERSION=1.20 \
    MINIMAP2_VERSION=2.30 \
    BUILD_DIR=/opt/build

RUN mkdir -p ${BUILD_DIR}/bin

# -------- Install HTSlib --------
RUN curl -fsSL -o htslib-${HTSLIB_VERSION}.tar.bz2 https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 \
    && tar -xjf htslib-${HTSLIB_VERSION}.tar.bz2 \
    && cd htslib-${HTSLIB_VERSION} \
    && ./configure --prefix=${BUILD_DIR} \
    && make -j"$(nproc)" && make install \
    && cd .. && rm -rf htslib-${HTSLIB_VERSION} htslib-${HTSLIB_VERSION}.tar.bz2

# -------- Install Samtools --------
RUN curl -fsSL -o samtools-${SAMTOOLS_VERSION}.tar.bz2 https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 \
    && tar -xjf samtools-${SAMTOOLS_VERSION}.tar.bz2 \
    && cd samtools-${SAMTOOLS_VERSION} \
    && ./configure --with-htslib=${BUILD_DIR} --prefix=${BUILD_DIR} \
    && make -j"$(nproc)" && make install \
    && cd .. && rm -rf samtools-${SAMTOOLS_VERSION} samtools-${SAMTOOLS_VERSION}.tar.bz2

# -------- Install Minimap2 --------
RUN curl -fsSL -o minimap2-${MINIMAP2_VERSION}.tar.bz2 https://github.com/lh3/minimap2/releases/download/v${MINIMAP2_VERSION}/minimap2-${MINIMAP2_VERSION}.tar.bz2 \
    && tar --no-same-owner -xjf minimap2-${MINIMAP2_VERSION}.tar.bz2 \
    && cd minimap2-${MINIMAP2_VERSION} \
    && arch="$(uname -m)"; \
       if [ "$arch" = "aarch64" ] || [ "$arch" = "arm64" ]; then \
           make -j"$(nproc)" aarch64=1; \
       else \
           make -j"$(nproc)"; \
       fi \
    && install -m 0755 minimap2 ${BUILD_DIR}/bin/minimap2 \
    && cd .. && rm -rf minimap2-${MINIMAP2_VERSION} minimap2-${MINIMAP2_VERSION}.tar.bz2

# -------- Add binaries to PATH --------
ENV PATH="${BUILD_DIR}/bin:${PATH}"

# -------- Create non-root user --------
RUN useradd -m -u 1000 -s /bin/bash user

# -------- Set working directory --------
WORKDIR /app
RUN chown user:user /app

# -------- Switch to non-root user --------
USER 1000

ENTRYPOINT ["/bin/bash"]
