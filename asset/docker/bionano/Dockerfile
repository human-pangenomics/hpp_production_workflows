FROM quay.io/masri2019/hpp_base:latest
MAINTAINER Mobin Asri, masri@ucsc.edu

RUN mkdir -p /home/apps/solve
RUN apt-get update && DEBIAN_FRONTEND="noninteractive" apt-get install -y python python-pip r-base
RUN pip install python-dateutil lxml pandas

## bionano solve
COPY ./tools/pipeline/1.0/HybridScaffold/1.0/scripts/fa2cmap_multi_color.pl /home/apps/solve
COPY ./tools/pipeline/1.0/RefAligner /home/apps/solve/RefAligner
COPY ./tools/pipeline/1.0/Pipeline /home/apps/solve/Pipeline

## UCSC convention is to work in /data
RUN mkdir -p /data
WORKDIR /data
