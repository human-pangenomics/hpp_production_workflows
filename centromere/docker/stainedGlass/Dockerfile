FROM continuumio/anaconda3
MAINTAINER Mobin Asri, masri@ucsc.edu

RUN mkdir -p /home/apps
RUN pip3 install scipy pandas matplotlib
RUN apt-get update
RUN apt-get install -y build-essential python3-dev autoconf automake make gcc perl zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev libssl-dev wget unzip curl libssl-dev libncurses5-dev libncursesw5-dev


# install mamba
RUN conda install -n base -c conda-forge mamba

# install snakemake
RUN mamba create -n snakemake -c conda-forge -c bioconda snakemake

# pull stainedGlass
RUN cd /home/apps && \
    git clone https://github.com/mrvollger/StainedGlass.git


# install samtools
RUN cd /home/apps && \
	wget https://github.com/samtools/samtools/releases/download/1.12/samtools-1.12.tar.bz2 && \
	tar -vxjf samtools-1.12.tar.bz2 && \
	rm -rf samtools-1.12.tar.bz2 && \
	cd samtools-1.12 && \
	make
ENV PATH="/home/apps/samtools-1.12:$PATH"

ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8

## UCSC convention is to work in /data
RUN mkdir -p /data
WORKDIR /data

#ENTRYPOINT ["conda", "run", "-n", "snakemake"]
