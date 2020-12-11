FROM tpesout/hpp_base:latest
MAINTAINER Trevor Pesout, tpesout@ucsc.edu

ENV PATH="/root/bin/python_3.6.0/bin/:${PATH}"
ENV PATH="/root/bin/samtools_1.9/:${PATH}"

RUN apt-get update && \
    apt-get install -y python-matplotlib

### QUAST
WORKDIR /opt
RUN git clone https://github.com/rlorigro/great_lengths.git && \
    cd great_lengths && \
    git checkout 796425fe8331f715bed762b591468a16402d4bae && \
    chmod +x great_lengths.py && \
    pip3 install matplotlib numpy && \
    ./great_lengths.py -h
ENV PATH="/opt/great_lengths/:${PATH}"

WORKDIR /data
