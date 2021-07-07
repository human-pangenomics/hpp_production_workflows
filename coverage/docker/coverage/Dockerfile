FROM quay.io/masri2019/hpp_base:latest
MAINTAINER Mobin Asri, masri@ucsc.edu

RUN mkdir -p /home/apps

COPY ./scripts /home/scripts
ENV DEPTH2COV_BIN=/home/scripts/bin/depth2cov
ENV CALC_MEAN_SD_PY=/home/scripts/calc_mean_sd.py
ENV COV2COUNTS_BIN=/home/scripts/bin/cov2counts
ENV FIT_MODEL_PY=/home/scripts/fit_model.py 
ENV PROJECT_BLOCKS_PY=/home/scripts/project_blocks.py
ENV ASSIGN_GAPS_PY=/home/scripts/assign_gaps.py
ENV FIT_MODEL_EXTRA_PY=/home/scripts/fit_model_extra.py
ENV FIND_BLOCKS_BIN=/home/scripts/bin/find_blocks_from_table

## UCSC convention is to work in /data
RUN mkdir -p /data
WORKDIR /data
