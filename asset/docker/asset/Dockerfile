FROM quay.io/masri2019/hpp_base:latest
MAINTAINER Mobin Asri, masri@ucsc.edu

RUN mkdir -p /home/apps

## asset
RUN cd /home/apps && git clone https://github.com/dfguan/asset && cd asset/src && make
ENV PATH="/home/apps/asset/src:${PATH}"

COPY ./scripts /home/scripts
ENV MERGE_BLOCKS_BASH=/home/scripts/merge_blocks.sh
ENV BLOCKS_STAT_BASH=/home/scripts/blocks_stat.sh
ENV BLOCKS_EXTENSION_V0_BASH=/home/scripts/blocks_extension_v0.sh
ENV BLOCKS_EXTENSION_V1_BASH=/home/scripts/blocks_extension_v1.sh

## UCSC convention is to work in /data
RUN mkdir -p /data
WORKDIR /data
