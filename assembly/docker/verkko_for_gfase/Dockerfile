
FROM humanpangenomics/verkko_hic@sha256:57fc0082db233941496c47c5dd61166169634e317af890500eafaab20366c2de
# FROM snakemake/snakemake:stable
MAINTAINER Julian Lucas, juklucas@ucsc.edu


RUN apt-get update && \
	apt-get install -y \
	build-essential \
	curl \
	automake \
	python3.8 \
	python3.8-dev \
	python3.8-distutils \
	python3-pip \
    autotools-dev \
    git \
    pigz \
    wget


## Install Rust (needed for Verkko)
RUN curl https://sh.rustup.rs -sSf | bash -s -- -y
ENV PATH="/root/.cargo/bin:${PATH}"
# RUN echo 'source $HOME/.cargo/env' >> $HOME/.bashrc

WORKDIR /opt

## Verkko
RUN git clone https://github.com/marbl/verkko.git \
	&& cd verkko/src \
	&& make -j 8

WORKDIR /data	