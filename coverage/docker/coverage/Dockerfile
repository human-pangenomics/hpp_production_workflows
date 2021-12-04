FROM quay.io/masri2019/hpp_base:latest
MAINTAINER Mobin Asri, masri@ucsc.edu

RUN mkdir -p /home/apps
RUN pip3 install scipy pandas matplotlib
RUN apt-get update
RUN apt-get install -y autoconf automake make gcc perl zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev libssl-dev wget unzip

#install Java
RUN cd /home && \
    wget https://download.oracle.com/java/17/latest/jdk-17_linux-x64_bin.tar.gz && \
    tar -xvzf jdk-17_linux-x64_bin.tar.gz
ENV PATH=$PATH:/home/jdk-17.0.1/bin

#install bcftools
RUN cd /home && \
        wget https://github.com/samtools/bcftools/releases/download/1.12/bcftools-1.12.tar.bz2 && \
        tar -vxjf bcftools-1.12.tar.bz2 && \
        rm -rf bcftools-1.12.tar.bz2 && \
        cd bcftools-1.12 && \
        make
ENV PATH=/home/bcftools-1.12:$PATH


#intstall IGV
RUN cd /home && \
    wget https://data.broadinstitute.org/igv/projects/downloads/2.11/IGV_2.11.3.zip &&\
    unzip IGV_2.11.3.zip
ENV PATH=$PATH:/home/IGV_2.11.3

#install hstlib  C API
RUN wget https://github.com/samtools/htslib/releases/download/1.13/htslib-1.13.tar.bz2 && \
    tar -xvjf htslib-1.13.tar.bz2 && \
    cd htslib-1.13 && \
    autoreconf -i && \
    ./configure && \
    make && \
    make install
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib


COPY ./scripts /home/scripts
ENV DEPTH2COV_BIN=/home/scripts/bin/depth2cov
ENV CALC_MEAN_SD_PY=/home/scripts/calc_mean_sd.py
ENV COV2COUNTS_BIN=/home/scripts/bin/cov2counts
ENV FIT_MODEL_PY=/home/scripts/fit_model.py
ENV PROJECT_BLOCKS_PY=/home/scripts/project_blocks.py
ENV ASSIGN_GAPS_PY=/home/scripts/assign_gaps.py
ENV FIT_MODEL_EXTRA_PY=/home/scripts/fit_model_extra.py
ENV FIND_BLOCKS_BIN=/home/scripts/bin/find_blocks_from_table
ENV SPLIT_CONTIGS_COV_BIN=/home/scripts/bin/split_contigs_cov 
ENV SPLIT_CONTIGS_COV_V2_BIN=/home/scripts/bin/split_contigs_cov_v2
ENV PDF_GENERATOR_PY=/home/scripts/pdf_generator.py
ENV SPLIT_BAM_BY_READNAME_BIN=/home/scripts/bin/split_bam_by_readname
ENV PHASE_READS_BIN=/home/scripts/bin/phase_reads
ENV PHASE_READS_DEBUG_BIN=/home/scripts/bin/phase_reads_debug
ENV CORRECT_BAM_BIN=/home/scripts/bin/correct_bam
ENV FILTER_ALT_READS_BIN=/home/scripts/bin/filter_alt_reads
ENV EXTRACT_READS_BIN=/home/scripts/bin/extract_reads
ENV INCREASE_MAPQ_BIN=/home/scripts/bin/increase_mapq


## UCSC convention is to work in /data
RUN mkdir -p /data
WORKDIR /data
