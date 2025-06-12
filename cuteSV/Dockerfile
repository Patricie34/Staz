FROM python:3.9-slim

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        git \
        build-essential \
    && rm -rf /var/lib/apt/lists/*

RUN addgroup cutesvgroup && \
    useradd -m -u 1000 -s /bin/bash -g cutesvgroup cutesvuser

WORKDIR /app
RUN chown -R cutesvuser:cutesvgroup /app

# Clone and install cuteSV before switching to non-root
RUN git clone https://github.com/tjiangHIT/cuteSV.git /app/cuteSV
WORKDIR /app/cuteSV
RUN python setup.py install

# Now drop to non-root user
USER 1000
