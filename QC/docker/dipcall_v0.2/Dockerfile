FROM tpesout/hpp_base:latest
MAINTAINER Julian Lucas, juklucas@ucsc.edu

### dipcall
WORKDIR /opt/dipcall
RUN wget https://github.com/lh3/dipcall/releases/download/v0.2/dipcall-0.2_x64-linux.tar.bz2 && \
    tar -jxf dipcall-0.2_x64-linux.tar.bz2 && \
    rm dipcall-0.2_x64-linux.tar.bz2

WORKDIR /data
