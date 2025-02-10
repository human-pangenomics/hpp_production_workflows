FROM quay.io/biocontainers/pandas:2.2.1
MAINTAINER Julian Lucas, juklucas@ucsc.edu

LABEL   \
    software="assign_chromosomes" \ 
    container="assign_chromosomes" \ 
    about.summary="Container to run assembly chromosome assignment script."
    

WORKDIR /opt/

RUN mkdir -p /opt/scripts
COPY ./scripts/chrom_alias_and_assignment.py /opt/scripts/chrom_alias_and_assignment.py
RUN chmod +x /opt/scripts/chrom_alias_and_assignment.py

ENV PATH="/opt/scripts:$PATH"

WORKDIR /data/