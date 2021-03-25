FROM tpesout/hpp_base:latest
MAINTAINER Trevor Pesout, tpesout@ucsc.edu

### blast ###
WORKDIR /opt/blast
RUN wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.11.0/ncbi-blast-2.11.0+-src.tar.gz && \
    tar -xvf ncbi-blast-2.11.0+-src.tar.gz && \
    rm ncbi-blast-2.11.0+-src.tar.gz && \
    cd ncbi-blast-2.11.0+-src/c++/ && \
    ./configure && \
    make -j 8
ENV PATH="$PATH:/root/bin/samtools_1.9:/opt/blast/ncbi-blast-2.11.0+-src/c++/ReleaseMT/bin"

### vecscreen ###
RUN apt-get install -y ncbi-tools-bin

### blast taxonomy db ###
WORKDIR /opt/blast/taxdb
RUN wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz && \
    tar xvf taxdb.tar.gz && \
    rm taxdb.tar.gz
ENV BLASTDB=/opt/blast/taxdb

WORKDIR /data
