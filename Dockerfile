# First dockerfile

FROM ubuntu:20.04

# Set the working directory in the container
WORKDIR /app

# Install the necessary packages without caching to keep the image size smaller
RUN pip install --no-cache-dir pandas scikit-allel matplotlib argparse matplotlib-venn
