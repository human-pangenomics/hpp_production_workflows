FROM mobinasri/bio_base:latest
MAINTAINER Mobin Asri, masri@ucsc.edu

RUN mkdir -p /home/apps
RUN pip3 install scipy pandas matplotlib
RUN apt-get update
RUN apt-get install -y build-essential python3-dev autoconf automake make gcc perl zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev libssl-dev wget unzip


#install veritymap
RUN apt install -y software-properties-common && \
    add-apt-repository -y ppa:ubuntu-toolchain-r/test && \
    apt install -y gcc-9 g++-9
ENV CC=/usr/bin/gcc-9
ENV CXX=/usr/bin/g++-9

RUN cd /home/apps && \
    wget https://github.com/ablab/VerityMap/archive/refs/tags/v2.1.1-alpha.zip && \
    unzip v2.1.1-alpha.zip && \
    cd VerityMap-2.1.1-alpha/ && \
    pip install -r requirements.txt && \
    cd veritymap && \
    make

ENV VERITY_MAP_PY="/home/apps/VerityMap-2.1.1-alpha/veritymap/main.py"

ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8

## UCSC convention is to work in /data
RUN mkdir -p /data
WORKDIR /data
