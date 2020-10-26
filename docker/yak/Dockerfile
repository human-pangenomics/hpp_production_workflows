FROM tpesout/hpp_base:latest
MAINTAINER Trevor Pesout, tpesout@ucsc.edu

WORKDIR /opt/yak
RUN wget https://github.com/lh3/yak/archive/v0.1.tar.gz && \
	tar xvf v0.1.tar.gz && \
	rm v0.1.tar.gz && \
	cd yak-0.1 && \
	make
ENV PATH="/opt/yak/yak-0.1/:$PATH"

WORKDIR /data
