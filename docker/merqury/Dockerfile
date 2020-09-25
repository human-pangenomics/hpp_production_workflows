FROM tpesout/hpp_base:latest
MAINTAINER Trevor Pesout, tpesout@ucsc.edu


RUN apt-get update && \
    apt-get -y install libc-dev bzip2 libcurl4-openssl-dev openjdk-8-jdk ant fort77 gfortran build-essential xorg-dev \
        liblzma-dev libblas-dev gcc-multilib gobjc++ aptitude libreadline-dev r-base

### R
# 3.2.1
WORKDIR /opt/R
RUN wget https://cran.r-project.org/src/base/R-3/R-3.6.1.tar.gz && \
    tar xvf R-3.6.1.tar.gz && \
    rm -r /opt/R/R-3.6.1.tar.gz && \
    cd R-3.6.1 && \
    ./configure --enable-utf8  --prefix /root/bin/R_3.6.1 && \
    make && \
    make install && \
  	/root/bin/R_3.6.1/bin/R -e "install.packages('argparse', repos ='http://cran.us.r-project.org', dependencies=TRUE)" \
  	/root/bin/R_3.6.1/bin/R -e "install.packages('ggplot2', repos ='http://cran.us.r-project.org', dependencies=TRUE)" \
  	/root/bin/R_3.6.1/bin/R -e "install.packages('scales', repos ='http://cran.us.r-project.org', dependencies=TRUE)"
ENV PATH="/root/bin/R_3.6.1/bin:${PATH}"

### bedtools
# 2.29.2
WORKDIR /opt/bedtools2
RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.29.2/bedtools-2.29.2.tar.gz && \
	tar xvf bedtools-2.29.2.tar.gz && \
	rm bedtools-2.29.2.tar.gz && \
	mv bedtools2 bedtools2_2.29.2 && \
	cd bedtools2_2.29.2 && \
	make && \
	mkdir -p /root/bin/bedtools_2.29.2 && \
	cp /opt/bedtools2/bedtools2_2.29.2/bin/* /root/bin/bedtools_2.29.2
ENV PATH="/root/bin/bedtools_2.29.2:${PATH}"

### IGVTools
# 2.3.98
WORKDIR /opt/igvtools
RUN wget http://data.broadinstitute.org/igv/projects/downloads/2.3/igvtools_2.3.98.zip && \
    unzip igvtools_2.3.98.zip && \
    rm igvtools_2.3.98.zip && \
    mv IGVTools IGVTools_2.3.98
ENV PATH="/opt/igvtools/IGVTools_2.3.98:${PATH}"

### Meryl
# 1.0
WORKDIR /opt/meryl/
RUN wget https://github.com/marbl/meryl/releases/download/v1.0/meryl-1.0.Linux-amd64.tar.xz && \
	tar xvf meryl-1.0.Linux-amd64.tar.xz && \
	rm meryl-1.0.Linux-amd64.tar.xz
ENV PATH="/opt/meryl/meryl-1.0/Linux-amd64/bin:${PATH}"


## Install OpenJDK-8
#RUN apt-get update && \
#    apt-get install -y  && \
#    apt-get clean;
#
## Fix certificate issues
#RUN apt-get update && \
#    apt-get install ca-certificates-java && \
#    apt-get clean && \
#    update-ca-certificates -f;

# merqury
WORKDIR /opt/merqury
RUN wget https://github.com/marbl/merqury/archive/v1.1.tar.gz && \
    tar xvf v1.1.tar.gz && \
    rm v1.1.tar.gz
ENV PATH="/opt/merqury/merqury-1.1:${PATH}"

# Final environment setup
ENV JAVA_HOME /usr/lib/jvm/java-8-openjdk-amd64/
ENV PATH="/root/bin/python_3.6.0/bin:${PATH}"
ENV MERQURY="/opt/merqury/merqury-1.1"

# issue generating pngs with Rscript in docker container without X11
RUN apt-get install -y xvfb gsfonts-x11 xfonts-base xfonts-scalable xfonts-100dpi xfonts-75dpi && \
    mv /root/bin/R_3.6.1/bin/Rscript /root/bin/R_3.6.1/bin/.Rscript && \
    echo '#!/bin/bash\nxvfb-run .Rscript $@' >/root/bin/R_3.6.1/bin/Rscript && \
    chmod +x /root/bin/R_3.6.1/bin/Rscript

WORKDIR /data