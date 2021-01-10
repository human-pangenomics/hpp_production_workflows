FROM quay.io/masri2019/hpp_base:latest
MAINTAINER Mobin Asri, masri@ucsc.edu

RUN apt-get update && apt-get install -y bwa

## UCSC convention is to work in /data
RUN mkdir -p /data
WORKDIR /data
