# ===============================
# 1️⃣ Stage 1: Build alleleCount
# ===============================
FROM ubuntu:20.04 AS builder

# Set environment variables
ENV OPT=/opt/wtsi-cgp \
    PATH=/opt/wtsi-cgp/bin:/opt/wtsi-cgp/biobambam2/bin:$PATH \
    PERL5LIB=/opt/wtsi-cgp/lib/perl5 \
    LD_LIBRARY_PATH=/opt/wtsi-cgp/lib \
    LC_ALL=en_US.UTF-8 \
    LANG=en_US.UTF-8 \
    VER_HTSLIB="1.11" \
    VER_LIBDEFLATE="v1.6"

# Install required system dependencies
RUN apt-get update && apt-get install -yq --no-install-recommends \
    build-essential apt-transport-https curl ca-certificates make \
    bzip2 gcc locales wget zlib1g-dev libbz2-dev liblzma-dev \
    libcurl4-gnutls-dev libncurses5-dev libtasn1-dev nettle-dev \
    libgmp-dev libp11-kit-dev && \
    locale-gen en_US.UTF-8 && \
    update-locale LANG=en_US.UTF-8 && \
    rm -rf /var/lib/apt/lists/*

# Build external dependencies
COPY build/opt-build.sh build/
RUN bash build/opt-build.sh $OPT

# Build alleleCount
COPY . .
RUN bash build/opt-build-local.sh $OPT

# ===============================
# 2️⃣ Stage 2: Create Minimal Image
# ===============================
FROM ubuntu:20.04

# Set environment variables
ENV OPT=/opt/wtsi-cgp \
    PATH=/opt/wtsi-cgp/bin:/opt/wtsi-cgp/biobambam2/bin:$PATH \
    PERL5LIB=/opt/wtsi-cgp/lib/perl5 \
    LD_LIBRARY_PATH=/opt/wtsi-cgp/lib \
    LC_ALL=en_US.UTF-8 \
    LANG=en_US.UTF-8

# Install runtime dependencies (minimal)
RUN apt-get update && apt-get install -yq --no-install-recommends \
    locales curl ca-certificates libperlio-gzip-perl bzip2 psmisc \
    time zlib1g liblzma5 libncurses5 p11-kit && \
    locale-gen en_US.UTF-8 && \
    update-locale LANG=en_US.UTF-8 && \
    rm -rf /var/lib/apt/lists/*

# Copy built alleleCount from the builder stage
COPY --from=builder $OPT $OPT

# Create a non-root user
RUN adduser --disabled-password --gecos '' ubuntu && chsh -s /bin/bash ubuntu

USER ubuntu
WORKDIR /home/ubuntu

CMD ["/bin/bash"]
