# First dockerfile

FROM python:3.8-slim

# Set the working directory in the container
WORKDIR /app

# Install the necessary packages without caching to keep the image size smaller
RUN pip install --no-cache-dir pandas scikit-allel matplotlib argparse matplotlib-venn
