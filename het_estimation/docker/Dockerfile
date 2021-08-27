FROM quay.io/masri2019/hpp_base:latest
MAINTAINER Mobin Asri, masri@ucsc.edu

RUN mkdir -p /home/apps

COPY ./scripts /home/scripts
ENV HET_ESTIMATOR_PY=/home/scripts/het_estimator.py

## UCSC convention is to work in /data
RUN mkdir -p /data
WORKDIR /data
