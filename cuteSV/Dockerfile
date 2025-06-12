FROM python:3.9-slim

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        git \
        build-essential \
        pip \
    && rm -rf /var/lib/apt/lists/*

RUN addgroup cutesvgroup && \
    useradd -m -u 1000 -s /bin/bash -g cutesvgroup cutesvuser

WORKDIR /app
RUN chown -R cutesvuser:cutesvgroup /app

# Clone the cuteSV repository
RUN git clone https://github.com/tjiangHIT/cuteSV.git /app/cuteSV

# Install using pip (modern method)
WORKDIR /app/cuteSV
RUN pip install .

# Drop to non-root user
USER 1000
