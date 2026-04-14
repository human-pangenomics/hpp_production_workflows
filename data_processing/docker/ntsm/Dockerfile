FROM ubuntu:22.04

MAINTAINER Julian Lucas, juklucas@ucsc.edu

LABEL base_image="ubuntu:22.04"
LABEL software="ntsm"
LABEL software.version="1.2.1"
LABEL about.home="https://github.com/JustinChu/ntsm"
LABEL about.license_file="https://github.com/JustinChu/ntsm/blob/master/LICENSE"
LABEL about.summary="ntsm identifies sample swaps and incongruent samples by counting k-mers in sequence data"

## Install build dependencies and runtime tools
RUN apt-get update \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y \
        git \
        gcc \
        g++ \
        make \
        autoconf \
        automake \
        libtool \
        zlib1g-dev \
        samtools \
        pigz \
    && rm -rf /var/lib/apt/lists/*

## Clone and build ntsm v1.2.1
RUN git clone --branch 1.2.1 --depth 1 https://github.com/JustinChu/ntsm.git /opt/ntsm_src \
    && cd /opt/ntsm_src \
    && ./autogen.sh \
    && ./configure --prefix=/opt/ntsm \
    && make install

## Copy pre-built reference data into install prefix
RUN cp -r /opt/ntsm_src/data /opt/ntsm/data

ENV PATH="/opt/ntsm/bin:${PATH}"

WORKDIR /data
