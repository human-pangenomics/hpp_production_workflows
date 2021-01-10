FROM quay.io/masri2019/hpp_base:latest
MAINTAINER Mobin Asri, masri@ucsc.edu

RUN apt-get update && apt-get install -y zip unzip curl
RUN mkdir -p /home/apps
## Winnowmap
RUN cd /home/apps && \
    wget https://github.com/marbl/Winnowmap/archive/v2.0.zip && \
    unzip v2.0.zip && \
    cd Winnowmap-2.0 && \
    make -j8
ENV PATH="/home/apps/Winnowmap-2.0/bin:${PATH}"

## minimap2
RUN cd /home/apps && \
    git clone https://github.com/lh3/minimap2 && \
    cd minimap2 && \
    make
ENV PATH="/home/apps/minimap2:${PATH}"
ENV PAFTOOLS_PATH="/home/apps/minimap2/misc/paftools.js"

## k8
RUN cd /home/apps && \
    curl -L https://github.com/attractivechaos/k8/releases/download/v0.2.4/k8-0.2.4.tar.bz2 | tar -jxf -
RUN cd /home/apps/k8-0.2.4/ && cp k8-`uname -s` k8
ENV PATH="/home/apps/k8-0.2.4:$PATH"

## UCSC convention is to work in /data
RUN mkdir -p /data
WORKDIR /data

