# ------------------------------------------------------------
#  Docker image: latest Minimap2 + latest stable Samtools
#  Ubuntu 22.04 base, non‑root user, /data as workdir
# ------------------------------------------------------------
FROM ubuntu:22.04

LABEL description="Docker image with the latest Minimap2 and Samtools (built from source)."

# -----------------------------------------------------------------
#  Environment
# -----------------------------------------------------------------
ENV DEBIAN_FRONTEND=noninteractive \
    # make sure the binaries we install go into the PATH
    PATH="/usr/local/bin:${PATH}"

# -----------------------------------------------------------------
#  Install system‑level build dependencies
# -----------------------------------------------------------------
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        ca-certificates \
        wget \
        curl \
        gnupg \
        git \
        build-essential \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        libcurl4-openssl-dev \
        libncurses5-dev libncursesw5-dev \
        libssl-dev \
        libgomp1 && \
    rm -rf /var/lib/apt/lists/*

# -----------------------------------------------------------------
#  Install Minimap2 (latest master branch)
# -----------------------------------------------------------------
WORKDIR /opt
RUN git clone --depth 1 https://github.com/lh3/minimap2.git && \
    cd minimap2 && \
    make -j$(nproc) && \
    cp minimap2 /usr/local/bin/ && \
    cd .. && rm -rf minimap2

# -----------------------------------------------------------------
#  Install Samtools (latest stable release) – built from source
# -----------------------------------------------------------------
# 1️⃣  Clone Samtools **with its htslib sub‑module**
# 2️⃣  Build & install
WORKDIR /opt
RUN git clone --depth 1 --branch master \
              --recurse-submodules https://github.com/samtools/samtools.git && \
    cd samtools && \
    ./configure --prefix=/usr/local && \
    make -j$(nproc) && \
    make install && \
    cd .. && rm -rf samtools

# -----------------------------------------------------------------
#  Create a non‑root user (UID 1000) and set ownership of /data
# -----------------------------------------------------------------
RUN groupadd -g 1000 app && \
    useradd -m -u 1000 -g app -s /bin/bash appuser

# -----------------------------------------------------------------
#  Working directory (where users will drop their files)
# -----------------------------------------------------------------
WORKDIR /data
RUN chown appuser:app /data

# -----------------------------------------------------------------
#  Switch to the non‑root user for the rest of the container life‑time
# -----------------------------------------------------------------
USER appuser
