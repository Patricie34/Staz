FROM ubuntu:22.04

# Set non-interactive mode
ENV DEBIAN_FRONTEND=noninteractive

# Install basic tools
RUN apt-get update && apt-get install -y \
    git \
    wget \
    curl \
    unzip \
    build-essential \
    python3 \
    python3-pip \
    python3-setuptools \
    bedtools \
    samtools \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Create python symlink
RUN ln -s /usr/bin/python3 /usr/bin/python

# ---------------------------
# Install non-B DNA GFA tool
# ---------------------------
RUN git clone https://github.com/abcsFrederick/non-B_gfa.git /opt/non-B_gfa && \
    cd /opt/non-B_gfa && \
    make && \
    cp gfa /usr/local/bin/ && \
    chmod +x /usr/local/bin/gfa

# ---------------------------
# Python utilities (optional)
# ---------------------------
RUN pip3 install --upgrade pip setuptools && \
    pip3 install pandas numpy pybedtools

# ---------------------------
# Create non-root user
# ---------------------------
RUN useradd -m -u 1000 -s /bin/bash user

# ---------------------------
# PATH exports
# ---------------------------
ENV PATH="/usr/local/bin:${PATH}"

# ---------------------------
# Default workdir and entrypoint
# ---------------------------
WORKDIR /data
RUN chown -R user:user /data
USER 1000
ENTRYPOINT ["/bin/bash"]
