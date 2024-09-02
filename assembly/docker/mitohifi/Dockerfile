FROM ghcr.io/marcelauliano/mitohifi@sha256:d8a72cac5add1d2836d9be4c2a50b32bb352699d1778d0ddcaaca79e1ce604ab
MAINTAINER Julian Lucas, juklucas@ucsc.edu

LABEL   \
    software="mitohifi" \ 
    container="mitohifi" \ 
    about.summary="Update released MitoHiFi container to Hifiasm 19.9"


# Hifiasm-0.19.9
RUN cd /opt/ && \
    wget https://github.com/chhylp123/hifiasm/archive/refs/tags/0.19.9.zip && \
    unzip 0.19.9.zip  && \
        cd hifiasm-0.19.9 && \ 
    make
ENV PATH="/opt/hifiasm-0.19.9:$PATH"

WORKDIR /tmp