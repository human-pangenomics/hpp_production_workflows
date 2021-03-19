FROM quay.io/masri2019/hpp_base:latest
MAINTAINER Mobin Asri, masri@ucsc.edu

RUN mkdir -p /home/apps

## asset
RUN cd /home/apps && git clone https://github.com/dfguan/asset && cd asset/src && make
ENV PATH="/home/apps/asset/src:${PATH}"

COPY ./scripts /home/scripts
ENV MERGE_BLOCKS_BASH=/home/scripts/merge_blocks.sh
ENV BLOCKS_STAT_V0_BASH=/home/scripts/blocks_stat_v0.sh
ENV BLOCKS_STAT_V1_BASH=/home/scripts/blocks_stat_v1.sh
ENV BLOCKS_STATS_PY=/home/scripts/blocks_stat.py
ENV BLOCKS_EXTENSION_V0_BASH=/home/scripts/blocks_extension_v0.sh
ENV BLOCKS_EXTENSION_V1_BASH=/home/scripts/blocks_extension_v1.sh
ENV DEPTH2COV_PY=/home/scripts/depth2cov.py
ENV CALC_MEAN_SD_PY=/home/scripts/calc_mean_sd.py
ENV COV2COUNTS_PY=/home/scripts/cov2counts.py
ENV COV_FREQ_DIST_PY=/home/scripts/coverage_frequency_distribution.py
ENV FIT_MODEL_PY=/home/scripts/fit_model.py 
ENV MAP_BLOCKS_PY=/home/scripts/map_blocks.py

## UCSC convention is to work in /data
RUN mkdir -p /data
WORKDIR /data
