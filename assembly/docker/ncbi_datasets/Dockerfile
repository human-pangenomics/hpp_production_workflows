FROM python:3.11-slim-buster
MAINTAINER Julian Lucas, juklucas@ucsc.edu

LABEL   \
    software="ncbi-datasets-cli" \ 
    container="ncbi-datasets-cli" \ 
    about.summary="ncbi datasets cli plus samtools/tabix/bcftools for compressing and indexing fasta files."


# Update package lists
RUN apt-get update && \
    apt-get install -y curl unzip samtools tabix bcftools && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* 


## install datasets & dataformat
WORKDIR /opt/
RUN curl -o datasets 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets' && \
    curl -o dataformat 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/dataformat' && \
    chmod +x datasets dataformat
ENV PATH=":/opt/:$PATH"

WORKDIR /data