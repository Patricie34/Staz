# First dockerfile

FROM ubuntu:20.04

# Set the working directory in the container
WORKDIR /app

# Install necessary packages without caching to keep the image size smaller
RUN apt-get update && apt-get install -y python3-pip && \
    pip3 install --no-cache-dir pandas scikit-allel matplotlib argparse matplotlib-venn plotly
