FROM python:3.11-slim-buster
MAINTAINER Julian Lucas, juklucas@ucsc.edu

LABEL \
    software="sequence_toolbox" \
    container="sequence_toolbox" \
    about.summary="Container with common bioinformatics tools including bedtools, minimap2, samtools, seqkit, and pigz"


# Update package lists
RUN apt-get update && \
    apt-get install -y wget git pigz bedtools samtools tabix bcftools && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* 

# Install general Python dependencies
RUN pip install --no-cache-dir --upgrade pip
RUN pip3 install --no-cache-dir numpy pysam pandas matplotlib seaborn biopython

## install minimap2
WORKDIR /opt/
RUN wget https://github.com/lh3/minimap2/releases/download/v2.28/minimap2-2.28_x64-linux.tar.bz2 && \
    tar -jxvf minimap2-2.28_x64-linux.tar.bz2
ENV PATH=":/opt/minimap2-2.28_x64-linux/:$PATH"

## install seqkit
WORKDIR /opt/seqkit
RUN wget https://github.com/shenwei356/seqkit/releases/download/v0.16.1/seqkit_linux_amd64.tar.gz && \
    tar -xvf seqkit_linux_amd64.tar.gz
ENV PATH=":/opt/seqkit/:$PATH"


WORKDIR /data