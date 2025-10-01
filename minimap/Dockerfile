# ------------------------------------------------------------
#  Dockerfile: samtools + minimap2 (latest) + non‑root user
# ------------------------------------------------------------
FROM ubuntu:22.04

LABEL maintainer="Your Name <you@example.com>"
LABEL description="Ubuntu 22.04 image with the latest Samtools & Minimap2 built from source, running as a non‑root user."

ENV DEBIAN_FRONTEND=noninteractive

# -----------------------------------------------------------------
#  Install build‑time system packages
# -----------------------------------------------------------------
# The back‑slashes below are *exactly* at the end of each line – no trailing spaces,
# no comments after them.
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        ca-certificates \
        wget \
        curl \
        gnupg \
        build-essential \                # gcc, make, etc.
        libncurses5-dev \
        libncursesw5-dev \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        libcurl4-gnutls-dev \
        libssl-dev \
        libgomp1 \
        git \                            # needed for cloning the sources
        python3 \
        python3-pip && \
    rm -rf /var/lib/apt/lists/*

# -----------------------------------------------------------------
#  Build Samtools (clone with sub‑modules)
# -----------------------------------------------------------------
ARG SAMTOOLS_REPO="https://github.com/samtools/samtools.git"
ARG SAMTOOLS_TAG="master"          # change to a tag (e.g. v1.21) for reproducibility

RUN git clone --depth 1 \
              --branch ${SAMTOOLS_TAG} \
              --recurse-submodules \
              ${SAMTOOLS_REPO} /opt/samtools && \
    cd /opt/samtools && \
    make -j$(nproc) && \
    make install && \
    cd / && \
    rm -rf /opt/samtools

# -----------------------------------------------------------------
#  Build Minimap2 (latest master – change TAG if you prefer a fixed version)
# -----------------------------------------------------------------
ARG MINIMAP2_REPO="https://github.com/lh3/minimap2.git"
ARG MINIMAP2_TAG="master"

RUN git clone --depth 1 \
              --branch ${MINIMAP2_TAG} \
              ${MINIMAP2_REPO} /opt/minimap2 && \
    cd /opt/minimap2 && \
    make -j$(nproc) && \
    cp minimap2 mm2-fast mm2-precise /usr/local/bin/ && \
    cd / && \
    rm -rf /opt/minimap2

# -----------------------------------------------------------------
#  Verify that the binaries are present (optional but helpful)
# -----------------------------------------------------------------
RUN samtools --version && minimap2 --version

# -----------------------------------------------------------------
#  Create a non‑root user
# -----------------------------------------------------------------
RUN groupadd -g 1000 app && \
    useradd -m -u 1000 -g app -s /bin/bash appuser

# -----------------------------------------------------------------
#  Working directory owned by the non‑root user
# -----------------------------------------------------------------
WORKDIR /app
RUN chown appuser:app /app

# -----------------------------------------------------------------
#  Switch to the non‑root user for the rest of the container life‑time
# -----------------------------------------------------------------
USER appuser

# -----------------------------------------------------------------
#  Default entrypoint – a shell so the image can be used as a toolbox
# -----------------------------------------------------------------
ENTRYPOINT ["/bin/bash"]
