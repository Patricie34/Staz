# ------------------------------------------------------------
#  Docker image: latest Minimap2 + latest stable Samtools
#  Ubuntu 22.04 base, non‑root user, /data as workdir
# ------------------------------------------------------------
FROM ubuntu:22.04

LABEL maintainer="your_email@example.com"
LABEL description="Docker image with the latest Minimap2 and Samtools (built from source)."

# -----------------------------------------------------------------
#  Environment
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
#  Install Samtools (latest stable release from the git repo)
# -----------------------------------------------------------------
# 1️⃣ Clone Samtools *with* its htslib sub‑module
# 2️⃣ Generate the configure script with autoreconf
# 3️⃣ Build & install
WORKDIR /opt
RUN git clone --depth 1 --branch master \
              --recurse-submodules https://github.com/samtools/samtools.git && \
    cd samtools && \
    autoreconf -i && \ 
    ./configure --prefix=/usr/local && \
    make -j$(nproc) && \
    make install && \
    cd .. && rm -rf samtools

# -----------------------------------------------------------------
#  Verify that the binaries are there (optional but helpful)
# -----------------------------------------------------------------
RUN samtools --version && minimap2 --version

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

# -----------------------------------------------------------------
#  Default entrypoint – a shell so the image can be used as a toolbox
# -----------------------------------------------------------------
ENTRYPOINT ["/bin/bash"]
