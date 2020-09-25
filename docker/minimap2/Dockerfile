FROM tpesout/hpp_base:latest
MAINTAINER Trevor Pesout, tpesout@ucsc.edu

WORKDIR /opt/minimap2
RUN wget https://github.com/lh3/minimap2/releases/download/v2.17/minimap2-2.17.tar.bz2 && \
	tar xvf minimap2-2.17.tar.bz2 && \
	rm minimap2-2.17.tar.bz2 && \
	cd minimap2-2.17 && \
	make
ENV PATH="/opt/minimap2/minimap2-2.17:$PATH"

WORKDIR /data
