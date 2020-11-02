FROM tpesout/hpp_base:latest
MAINTAINER Trevor Pesout, tpesout@ucsc.edu

ENV MERYL_GIT_COMMIT="1fdad3e85652f033d4abea9045ff164634b645a8"


RUN apt-get update && \
    apt-get -y install libc-dev bzip2 libcurl4-openssl-dev openjdk-8-jdk ant fort77 gfortran build-essential xorg-dev \
        liblzma-dev libblas-dev gcc-multilib gobjc++ aptitude libreadline-dev r-base xvfb gsfonts-x11 xfonts-base \
        xfonts-scalable xfonts-100dpi xfonts-75dpi

### R
# 3.2.1
WORKDIR /opt/R
RUN wget https://cran.r-project.org/src/base/R-3/R-3.6.3.tar.gz && \
    tar xvf R-3.6.3.tar.gz && \
    rm -r /opt/R/R-3.6.3.tar.gz && \
    cd R-3.6.3 && \
    ./configure --enable-utf8  --prefix /root/bin/R_3.6.3 && \
    make && \
    make install && \
  	/root/bin/R_3.6.3/bin/R -e "install.packages('argparse', repos ='http://cran.us.r-project.org', dependencies=TRUE)" \
  	/root/bin/R_3.6.3/bin/R -e "install.packages('ggplot2', repos ='http://cran.us.r-project.org', dependencies=TRUE)" \
  	/root/bin/R_3.6.3/bin/R -e "install.packages('scales', repos ='http://cran.us.r-project.org', dependencies=TRUE)"
ENV PATH="/root/bin/R_3.6.3/bin:${PATH}"

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
RUN git clone https://github.com/marbl/meryl.git && \
    cd meryl && \
    git fetch && \
    git checkout $MERYL_GIT_COMMIT && \
    cd src && \
    make
ENV PATH="/opt/meryl/meryl/build/bin:${PATH}"


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
ENV MERQURY="/opt/merqury/merqury-1.1"
ENV PATH="/root/bin/python_3.6.0/bin:${PATH}"
ENV PATH="/root/bin/samtools_1.9:${PATH}"

# issue generating pngs with Rscript in docker container without X11
#RUN mv /root/bin/R_3.6.1/bin/Rscript /root/bin/R_3.6.1/bin/.Rscript && \
#    echo '#!/bin/bash\nxvfb-run .Rscript $@' >/root/bin/R_3.6.1/bin/Rscript && \
#    chmod +x /root/bin/R_3.6.1/bin/Rscript

WORKDIR /data