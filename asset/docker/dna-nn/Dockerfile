FROM quay.io/masri2019/hpp_base:latest
MAINTAINER Mobin Asri, masri@ucsc.edu

RUN apt-get update
RUN mkdir -p /home/apps
RUN cd /home/apps/ && \ 
	git clone https://github.com/lh3/dna-nn && \
	cd dna-nn && \
	make

ENV PATH="/home/apps/dna-nn:$PATH"

## UCSC convention is to work in /data
RUN mkdir -p /data
WORKDIR /data
