FROM continuumio/miniconda3:4.10.3

# OS requirements
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends --no-install-suggests -y \
    build-essential && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Set the conda channels
RUN conda config --add channels bioconda && \
    conda config --add channels conda-forge

# Add dependencies
RUN conda install -y \
    python==3.8 \
    rust==1.60.0 \
    dendropy==4.5.2 \
    numpy==1.21.1 && \
    conda clean -afy

RUN python -m pip install \
    setuptools_rust==1.2.0 \
    setuptools \
    wheel

# Add the scripts
RUN mkdir -p /phylodm
COPY . /phylodm
WORKDIR /phylodm
RUN python -m pip install .

# Run unit tests
RUN mkdir -p /test
COPY ./test /test
WORKDIR /test
RUN python -m unittest discover test
