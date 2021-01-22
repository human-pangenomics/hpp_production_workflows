FROM quay.io/masri2019/hpp_base:latest
MAINTAINER Mobin Asri, masri@ucsc.edu

RUN mkdir -p /home/apps/pstools

## asset
RUN cd /home/apps/pstools && wget https://pstools.s3.us-east-2.amazonaws.com/pstools && chmod 777 pstools
ENV PATH="/home/apps/pstools:${PATH}"

## UCSC convention is to work in /data
RUN mkdir -p /data
WORKDIR /data
