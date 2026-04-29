FROM python:3.11-slim

# System dependencies
RUN apt-get update && apt-get install -y \
    gcc \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    zlib1g-dev \
    procps \
    && rm -rf /var/lib/apt/lists/*

# Create non-root user
RUN groupadd -g 1000 user && \
    useradd -m -u 1000 -g 1000 -s /bin/bash user

# Set working directory
WORKDIR /app

# Copy requirements first
COPY requirements.txt .

# Install Python deps (still as root, system-wide)
RUN pip install --no-cache-dir -r requirements.txt

# Copy project files
COPY . .

# Fix permissions so user can write outputs
RUN chown -R user:user /app

# Switch to non-root user
USER user

# Default command
CMD ["python", "your_script.py"]
