FROM ubuntu:20.04
MAINTAINER Julian Lucas, juklucas@ucsc.edu

LABEL    software="verkko_hic" \ 
    container="verkko-hic" \ 
    about.summary="HiC Phasing For Verkko Assembly Runs" \ 
    about.home="https://github.com/Dmitry-Antipov/verkkohic"


RUN apt-get update && apt-get -y upgrade \
  && apt-get install -y --no-install-recommends \
    git \
    wget \
    g++ \
    pigz \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*


## Install Verkko, Mashmap, & Networkx w/ Miniconda
ENV PATH="/root/miniconda3/bin:${PATH}"
ARG PATH="/root/miniconda3/bin:${PATH}"
ADD environment.yml /tmp/environment.yml

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && mkdir /root/.conda \
    && bash Miniconda3-latest-Linux-x86_64.sh -b \
    && rm -f Miniconda3-latest-Linux-x86_64.sh \
    && conda init bash \
    && . /root/.bashrc \
    && conda update conda \
    && conda env create -f /tmp/environment.yml


## pstools
WORKDIR /opt
RUN wget https://pstools.s3.us-east-2.amazonaws.com/pstools_1 \
    && chmod +x pstools_1
ENV PATH="/opt:${PATH}"


## verkko hic
RUN git clone https://github.com/Dmitry-Antipov/verkkohic.git && \
    cd verkkohic && \
    git reset --hard c2669d7df7de3af09853727fb1fcbf650cedd02c

ENV PATH="/opt/verkkohic:${PATH}"

## Verkko HiC requires the following parameters to be set:
ENV PSTOOLS="/opt"
ENV VERKKO="/root/miniconda3/envs/verkko_hic"


## Make sure to enter the environment when starting
RUN echo "source activate $(head -1 /tmp/environment.yml | cut -d' ' -f2)" > ~/.bashrc
ENV PATH /opt/conda/envs/$(head -1 /tmp/environment.yml | cut -d' ' -f2)/bin:$PATH

WORKDIR /data
ENTRYPOINT [ "/bin/bash", "-l", "-c" ]

