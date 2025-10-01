# ============================================================
#  Minimal Genomics Toolbox – HTSlib + Samtools + Minimap2
#  (ncurses development packages added so Samtools builds)
# ============================================================

# ------------------------------------------------------------
#  Stage 1 – compile everything
# ------------------------------------------------------------
ARG UBUNTU_VERSION=24.04          # change to 22.04 if you prefer
ARG HTSLIB_VERSION=1.21           # adjust to the release you need
ARG SAMTOOLS_VERSION=1.21         # must match the HTSlib above
ARG MINIMAP2_VERSION=2.30         # latest at time of writing

FROM ubuntu:${UBUNTU_VERSION} AS builder

# ----- re‑declare the version arguments for this stage ----------
ARG HTSLIB_VERSION
ARG SAMTOOLS_VERSION
ARG MINIMAP2_VERSION

# ----- install only the packages needed for building ----------
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        ca-certificates \
        wget \
        curl \
        git \
        build-essential \          # gcc, g++, make, …
        autoconf \
        automake \
        libtool \
        pkg-config \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        libcurl4-openssl-dev \
        libssl-dev \
        # ***  <-- added ncurses development libs ***
        libncurses5-dev \
        libncursesw5-dev && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /opt/build

# ------------------------------------------------------------
#  HTSlib (foundation for Samtools)
# ------------------------------------------------------------
RUN wget -O htslib-${HTSLIB_VERSION}.tar.bz2 \
        https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 && \
    tar -xjf htslib-${HTSLIB_VERSION}.tar.bz2 && \
    cd htslib-${HTSLIB_VERSION} && \
    ./configure --prefix=/opt/build && \
    make -j$(nproc) && \
    make install && \
    cd /opt/build && rm -rf htslib-${HTSLIB_VERSION}*

# ------------------------------------------------------------
#  Samtools (built against the HTSlib we just installed)
# ------------------------------------------------------------
RUN wget -O samtools-${SAMTOOLS_VERSION}.tar.bz2 \
        https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    tar -xjf samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    cd samtools-${SAMTOOLS_VERSION} && \
    ./configure --prefix=/opt/build --with-htslib=/opt/build && \
    make -j$(nproc) && \
    make install && \
    cd /opt/build && rm -rf samtools-${SAMTOOLS_VERSION}*

# ------------------------------------------------------------
#  Minimap2 (independent binary)
# ------------------------------------------------------------
RUN wget -O minimap2-${MINIMAP2_VERSION}.tar.bz2 \
        https://github.com/lh3/minimap2/releases/download/v${MINIMAP2_VERSION}/minimap2-${MINIMAP2_VERSION}.tar.bz2 && \
    tar -xjf minimap2-${MINIMAP2_VERSION}.tar.bz2 && \
    cd minimap2-${MINIMAP2_VERSION} && \
    make -j$(nproc) && \
    install -m 0755 minimap2 /opt/build/bin/minimap2 && \
    cd /opt/build && rm -rf minimap2-${MINIMAP2_VERSION}*

# Refresh the dynamic‑linker cache so the binaries can find libhts.so
RUN ldconfig

# ------------------------------------------------------------
#  Stage 2 – tiny runtime image
# ------------------------------------------------------------
FROM ubuntu:${UBUNTU_VERSION}

ENV DEBIAN_FRONTEND=noninteractive \
    LANG=C.UTF-8 \
    LC_ALL=C.UTF-8

# ---- runtime libraries required by the compiled tools ----
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        ca-certificates \
        libcurl4 \
        libssl3 \
        libzstd1 \
        libbz2-1.0 \
        liblzma5 \
        libdeflate0 \
        libncursesw6 && \
    rm -rf /var/lib/apt/lists/*

# ---- copy the compiled tools and libraries from the builder stage ----
COPY --from=builder /opt/build /usr/local

# Ensure the dynamic linker sees the libraries we just copied
RUN ldconfig

# ---- create a non‑root user ---------------------------------
ARG USER=worker
ARG UID=2000
ARG GID=2000
RUN groupadd -g ${GID} ${USER} && \
    useradd -m -u ${UID} -g ${GID} -s /bin/bash ${USER}

# ---- environment -------------------------------------------
ENV PATH="/usr/local/bin:${PATH}"
WORKDIR /data
USER ${USER}

# ---- default command (interactive shell) --------------------
CMD ["bash"]
