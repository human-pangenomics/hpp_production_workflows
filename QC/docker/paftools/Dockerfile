FROM tpesout/hpp_base:latest
MAINTAINER Julian Lucas, juklucas@ucsc.edu

WORKDIR /opt

RUN git clone https://github.com/lh3/minimap2 && \
	cd minimap2 && \
	git reset --hard 8170693de39b667d11c8931d343c94a23a7690d2 && \
	make

ENV PATH="/opt/minimap2:/opt/minimap2/misc:$PATH"


WORKDIR /opt/k8
RUN wget https://github.com/attractivechaos/k8/releases/download/v0.2.4/k8-0.2.4.tar.bz2 && \
    tar xvf k8-0.2.4.tar.bz2 && \
    rm k8-0.2.4.tar.bz2 && \
    cp k8-0.2.4/k8-`uname -s` k8-0.2.4/k8
ENV PATH="/opt/k8/k8-0.2.4:$PATH"


WORKDIR /data
