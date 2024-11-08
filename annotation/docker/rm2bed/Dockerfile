FROM quay.io/biocontainers/bioframe:0.3.3--pyhdfd78af_0
MAINTAINER Julian Lucas, juklucas@ucsc.edu

LABEL   \
    software="RM2BED" \ 
    container="rm2bed" \ 
    about.summary="Container to run RepeatMasker's RM2BED.py script."


WORKDIR /opt/
RUN wget https://raw.githubusercontent.com/Dfam-consortium/RepeatMasker/c08c50dfa026a2422150054033878f8e729ca5ea/util/RM2Bed.py && \
    chmod +x RM2Bed.py

ENV PATH=":/opt/:$PATH"    

WORKDIR /data/