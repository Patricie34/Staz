FROM continuumio/miniconda3

ENV DEBIAN_FRONTEND=noninteractive

# Create group and user with gid and uid 1000
RUN groupadd -g 1000 nanovargroup && \
    useradd -m -u 1000 -g nanovargroup -s /bin/bash nanovar

# Set working directory to /app/nanovar
WORKDIR /app/nanovar

# Install bioinformatics tools + python in conda env
RUN conda create -n nanovar-env -c conda-forge -y \
    python=3.11 \
    bedtools=2.26.0 \
    samtools=1.3.0 \
    minimap2=2.17 \
    && conda clean -afy

# Add requirements.txt inline with pinned Python libs
RUN echo "absl-py==2.1.0
astunparse==1.6.3
beautifulsoup4==4.12.3
biopython==1.84
bs4==0.0.2
cachetools==5.5.0
certifi==2024.12.14
charset-normalizer==3.4.1
contourpy==1.3.1
cycler==0.12.1
flatbuffers==24.12.23
fonttools==4.55.3
gast==0.6.0
google-auth==2.37.0
google-auth-oauthlib==1.2.1
google-pasta==0.2.0
grpcio==1.69.0
h5py==3.12.1
idna==3.10
keras==2.15.0
kiwisolver==1.4.8
libclang==18.1.1
Markdown==3.7
MarkupSafe==3.0.2
matplotlib==3.10.0
ml-dtypes==0.3.2
nanovar==1.8.3
natsort==8.4.0
numpy==1.26.4
oauthlib==3.2.2
opt_einsum==3.4.0
packaging==24.2
pandas==2.2.3
pillow==11.1.0
pip==24.3.1
protobuf==4.25.5
pyasn1==0.6.1
pyasn1_modules==0.4.1
pybedtools==0.11.0
pyparsing==3.2.1
pysam==0.22.1
python-dateutil==2.9.0.post0
pytz==2024.2
requests==2.32.3
requests-oauthlib==2.0.0
rsa==4.9
scipy==1.15.0
setuptools==75.6.0
six==1.17.0
soupsieve==2.6
tensorboard==2.15.2
tensorboard-data-server==0.7.2
tensorflow-cpu==2.15.1
tensorflow-estimator==2.15.0
tensorflow-io-gcs-filesystem==0.37.1
termcolor==2.5.0
typing_extensions==4.12.2
tzdata==2024.2
urllib3==2.3.0
Werkzeug==3.1.3
wheel==0.45.1
wrapt==1.14.1
" > requirements.txt

# Install Python dependencies inside the conda env
RUN /bin/bash -c "source activate nanovar-env && pip install --no-cache-dir -r requirements.txt"

# Set conda environment PATH
ENV PATH /opt/conda/envs/nanovar-env/bin:$PATH

# Change ownership of /app/nanovar to user/group
RUN chown -R nanovar:nanovargroup /app/nanovar

# Switch to non-root user
USER 1000

WORKDIR /app/nanovar
