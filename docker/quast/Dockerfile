FROM tpesout/hpp_base:latest
MAINTAINER Trevor Pesout, tpesout@ucsc.edu

ENV PATH="/root/bin/python_2.7.13/bin:${PATH}"

RUN apt-get update && \
    apt-get install -y libfreetype6-dev libpng-dev python-matplotlib

# pkg-config python-setuptools libsasl2-dev libssl-dev


### QUAST
WORKDIR /root/tools/quast
RUN wget https://downloads.sourceforge.net/project/quast/quast-5.0.2.tar.gz && \
    tar -xzf quast-5.0.2.tar.gz && \
    rm quast-5.0.2.tar.gz && \
    cd quast-5.0.2 && \
    pip install joblib==0.14.1 && \
    python setup.py install


WORKDIR /data
