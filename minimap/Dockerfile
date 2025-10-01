# ------------------------------------------------------------
#  Docker image: latest Minimap2 + latest HTSlib + latest Samtools
#  Ubuntu 22.04 base, non‑root user, /data as workdir
# ------------------------------------------------------------
FROM ubuntu:22.04

LABEL description="Docker image with Minimap2, HTSlib and Samtools (all built from source)."

# -----------------------------------------------------------------
#  Global environment
# -----------------------------------------------------------------
ENV DEBIAN_FRONTEND=noninteractive \
    PATH="/usr/local/bin:${PATH}"

# -----------------------------------------------------------------
#  Install system‑level build dependencies (including autotools)
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
        libgomp1 \
        autoconf \
        automake \
        libtool \
        pkg-config && \
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
#  1️⃣ Build and install HTSlib (latest master)
# -----------------------------------------------------------------
WORKDIR /opt
RUN git clone --depth 1 --branch master \
              --recurse-submodules https://github.com/samtools/htslib.git && \
    cd htslib && \
    autoreconf -i && \
    ./configure --prefix=/usr/local && \
    make -j$(nproc) && \
    make install && \
    cd .. && rm -rf htslib

# -----------------------------------------------------------------
#  2️⃣ Build and install Samtools (latest master) – link to the HTSlib we just built
# -----------------------------------------------------------------
WORKDIR /opt
RUN git clone --depth 1 --branch master \
              --recurse-submodules https://github.com/samtools/samtools.git && \
    cd samtools && \
    # Ensure the submodule is fully present (no‑op if already complete)
    git submodule update --init --recursive && \
    autoreconf -i && \
    ./configure --prefix=/usr/local --with-htslib=/usr/local && \
    make -j$(nproc) && \
    make install && \
    cd .. && rm -rf samtools

# -----------------------------------------------------------------
#  Verify that the three binaries are available (optional but helpful)
# -----------------------------------------------------------------
RUN samtools --version && \
    minimap2 --version && \
    htslib --version

# -----------------------------------------------------------------
#  Create a non‑root user (UID 1000) and give it ownership of /data
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
