FROM tpesout/hpp_base:latest
MAINTAINER Trevor Pesout, tpesout@ucsc.edu

ENV PATH="/root/bin/python_3.6.0/bin/:${PATH}"

RUN apt-get update && \
    apt-get install -y pkg-config libfreetype6-dev libpng-dev python-matplotlib python-simplejson


### QUAST
WORKDIR /root/tools/quast
RUN wget https://downloads.sourceforge.net/project/quast/quast-5.0.2.tar.gz && \
    tar -xzf quast-5.0.2.tar.gz && \
    rm quast-5.0.2.tar.gz && \
    cd quast-5.0.2 && \
    pip3 install --upgrade setuptools matplotlib simplejson joblib==0.14.1 && \
    ./setup.py install && \
    ./quast.py -h


WORKDIR /data
