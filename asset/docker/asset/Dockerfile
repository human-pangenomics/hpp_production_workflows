FROM quay.io/masri2019/hpp_base:latest
MAINTAINER Mobin Asri, masri@ucsc.edu

RUN pip3 install numpy
RUN mkdir -p /home/apps

## asset
RUN cd /home/apps && git clone https://github.com/dfguan/asset && cd asset/src && make
ENV PATH="/home/apps/asset/src:${PATH}"

COPY ./scripts /home/scripts
ENV MERGE_BLOCKS_BASH=/home/scripts/merge_blocks.sh
ENV MERGE_BLOCKS_BASH_2=/home/scripts/merge_blocks_2.sh
ENV BLOCKS_STAT_V0_BASH=/home/scripts/blocks_stat_v0.sh
ENV BLOCKS_STAT_V1_BASH=/home/scripts/blocks_stat_v1.sh
ENV BLOCKS_STAT_V2_BASH=/home/scripts/blocks_stat_v2.sh
ENV BLOCKS_STAT_PY=/home/scripts/blocks_stat.py
ENV BLOCKS_EXTENSION_V0_BASH=/home/scripts/blocks_extension_v0.sh
ENV BLOCKS_EXTENSION_V1_BASH=/home/scripts/blocks_extension_v1.sh

## UCSC convention is to work in /data
RUN mkdir -p /data
WORKDIR /data
