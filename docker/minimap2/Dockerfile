FROM tpesout/hpp_base:latest
MAINTAINER Trevor Pesout, tpesout@ucsc.edu

WORKDIR /opt/minimap2
RUN wget https://github.com/lh3/minimap2/releases/download/v2.17/minimap2-2.17.tar.bz2 && \
	tar xvf minimap2-2.17.tar.bz2 && \
	rm minimap2-2.17.tar.bz2 && \
	cd minimap2-2.17 && \
	make
ENV PATH="/opt/minimap2/minimap2-2.17:/opt/minimap2/minimap2-2.17/misc:$PATH"

WORKDIR /opt/k8
RUN wget https://github.com/attractivechaos/k8/releases/download/v0.2.4/k8-0.2.4.tar.bz2 && \
    tar xvf k8-0.2.4.tar.bz2 && \
    rm k8-0.2.4.tar.bz2 && \
    cp k8-0.2.4/k8-`uname -s` k8-0.2.4/k8
ENV PATH="/opt/k8/k8-0.2.4:$PATH"

WORKDIR /data
