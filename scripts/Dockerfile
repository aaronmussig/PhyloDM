FROM continuumio/miniconda3:4.10.3

# OS requirements
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends --no-install-suggests -y \
    time \
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
    seaborn==0.11.2 \
    matplotlib-base==3.5.1 \
    pandas==1.4.2 \
    tqdm==4.64.0 \
    numpy==1.21.1 && \
    conda clean -afy

RUN python -m pip install \
    setuptools_rust==1.2.0 \
    setuptools \
    wheel

# Create a bind point for RDM
RUN mkdir -p /io

# Add the scripts
RUN mkdir -p /phylodm
COPY . /phylodm
WORKDIR /phylodm
RUN python -m pip install .

# Set the entrypoint
ENTRYPOINT ["python", "/phylodm/scripts/phylodm_perf.py"]
