FROM ubuntu:20.04
MAINTAINER Trevor Pesout, tpesout@ucsc.edu


RUN apt-get update && \
    DEBIAN_FRONTEND="noninteractive" apt-get install -y vim git wget build-essential python3 python3-dev python3-pip python3-matplotlib tabix

RUN pip3 install numpy matplotlib cython==0.29

WORKDIR /opt/
RUN git clone https://github.com/whatshap/whatshap.git && \
    cd whatshap && \
    git checkout e56c579d18d628270bc96b9411f302f16ed7758d && \
    python3 setup.py install

COPY chr_lengths /root/chr_lengths
COPY merge_results.py /root/merge_results.py

WORKDIR /data
